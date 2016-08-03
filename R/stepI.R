## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-04-08)

stepI <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  if ( x@locus@kind == "undefined" ) stop("undefined locus not allowed")
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- x@taxon@tip.rank
  msa.tab <- paste(tip.rank, gsub("^_", "", gene), sep = "_")
  block.max.dist <- x@params@block.max.dist
  max.bp <- x@params@max.bp * 1.5
  min.n.seq <- x@params@min.n.seq
  gblocks.exe <- x@mask.exe
  b1 <- x@params@gb1
  b2 <- x@params@gb2
  b3 <- x@params@gb3
  b4 <- x@params@gb4
  b5 <- x@params@gb5
  
  ## iniate logfile
  ## --------------
  logfile <- paste(gene, "stepI.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP I: detecting nucleotide positions of uncertain homology\n",
       paste("\n.. locus:", x@locus@sql),
       file = logfile)
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
  
  ## check if stepC has been run
  ## ---------------------------
  status <- paste("SELECT DISTINCT status",
                  "FROM", acc.tab)
  status <- dbGetQuery(conn, status)
  if ( "raw" %in% status$status ){
    dbDisconnect(conn)
    stop("stepC has not been run")
  }
  
  ## check if msa table exists
  ## -------------------------
  if ( !dbExistsTable(conn, msa.tab) ){
    dbDisconnect(conn)
    slog("\nWARNING: table", msa.tab, "does not exist!\n", file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP I finished after", round(td, 2), attr(td, "units"), "\n",
         file = logfile)
    return()
  }
  
  ## read alignment
  ## --------------
  slog("\n.. reading alignment", file = logfile)
  #   if ( x@params@parallel ){
  #     spec <- dbGetQuery(conn, paste("SELECT", tip.rank, "FROM", msa.tab))[, tip.rank]
  #     spec <- paste("^", spec, "$", sep = "")
  #     n <- length(spec)
  #     id <- seq(from = 1, to = n, by = ceiling(n/x@params@cpus))
  #     id <- data.frame(from = id, to = c(id[-1] - 1, n))
  #     spec <- apply(id, 1, function(i, a) paste(a[i[1]:i[2]], collapse = "|"), a = spec)
  #     sfInit(parallel = TRUE, cpus = x@params@cpus, 
  #            type = x@params@cluster.type)
  #     sfLibrary("megaptera", character.only = TRUE)
  #     sfLibrary("seqinr", character.only = TRUE)
  #     megProj <- x
  #     sfExport("spec", "msa.tab", "megProj")
  #     a <- sfLapply(x = spec, fun = dbPReadDNA, conn = megProj, tab.name = msa.tab, regex = TRUE)
  #     a <- do.call(rbind, a)
  #   } else {
  a <- dbReadDNA(x, msa.tab, taxon = ".+", regex = TRUE,
                 ignore.excluded = TRUE, blocks = "split")
  # }
  if ( is.null(a) ) {
    dbDisconnect(conn)
    slog("\nWARNING: no sequences conform to current parameter setting\n", file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP I finished after", round(td, 2), attr(td, "units"), "\n",
         file = logfile)
    return()
  }
  if ( is.matrix(a) ) a <- list(a)
  if ( !all(sapply(a, is.matrix)) ){
    stop("stepG has not been called yet")
  }
  slog("\n.. number of blocks:", length(a), file = logfile)
  n <- sapply(a, nrow)
  slog("\n.. number of taxa in blocks:", 
       paste(n, collapse = " - "), file = logfile)
  #   if ( nrow(a) == 1 ) {
  #     dbDisconnect(conn)
  #     slog("\n", file = logfile)
  #     dbWriteMSA(x, dna = a, masked = TRUE)
  #     return()
  #   }
  
  ## write files
  ## -----------
  #   slog("\n.. write alignment to files ..", file = logfile)
  #   write.phy(a, paste(gene, "phy", sep = "."))
  #   rownames(a) <- gsub("-", "_", rownames(a))
  #   write.nex(a, paste(gene, "nex", sep = "."))
  
  ## masking of poorly aligned nucleotides
  ## -----------------------------------
  slog("\n.. masking poorly aligned regions ..\n", file = logfile)
  a[n > 1] <- lapply(a[n > 1], gblocks, exec = gblocks.exe,
                     b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5) # with least conservative default
  a <- lapply(a, deleteEmptyCells, quiet = TRUE)
  
  ## sort alignment taxonomically
  ## ----------------------------
  if ( inherits(x@taxon, "taxonGuidetree") ){
    gt <- comprehensiveGuidetree(x, tip.rank = tip.rank, subset = a)
  } else { 
    gt <- dbReadTaxonomy(x, subset = a)
    gt <- tax2tree(gt, tip.rank = tip.rank)
  }
  matchAlignment <- function(ali, phy){
    if ( nrow(ali) > 2 ){
      nt <- setdiff(phy$tip.label, rownames(ali))
      phy <- ladderize(drop.tip(phy, nt))
      ali[match(phy$tip.label, rownames(ali)), ]
    }
    ali
  }
  a <- lapply(a, matchAlignment, phy = gt)
  
  ## write to database -- serially (or in parallel)
  ## ----------------------------------------------
  slog("\n.. write alignment/blocks to database ..", file = logfile)
  lapply(a, dbWriteMSA, megapteraProj = x, masked = TRUE)
  #   if ( x@params@parallel ){
  #     if ( nrow(a) > 2000 ){
  #       id <- seq(from = 1, to = nrow(a), by = ceiling(nrow(a)/x@params@cpus))
  #       id <- data.frame(from = id, to = c(id[-1] - 1, nrow(a)))
  #       aa <- apply(id, 1, function(i, a) a[i[1]:i[2], ], a = a)
  #       sfLibrary("seqinr", character.only = TRUE)
  #       sfExport("aa")
  #       sfLapply(x = aa, fun = dbWriteMSA, 
  #                megapteraProj = megProj, masked = TRUE)
  #     } else {
  #       dbWriteMSA(x, dna = a, masked = TRUE)
  #     }
  #     sfStop()
  #   } else {
  #     dbWriteMSA(x, dna = a, masked = TRUE)
  #   }
  
  ## write files
  ## -----------
  slog("\n.. write alignment to files ..", file = logfile)
  if ( length(a) == 1 ) {
    write.phy(a[[1]], paste(gene, "masked.phy", sep = "-"))
    rownames(a[[1]]) <- gsub("-", "_", rownames(a[[1]]))
    write.nex(a[[1]], paste(gene, "masked.nex", sep = "-"))
  } else {
    for ( i in 1:length(a) ){
      b <- paste("block", i, sep = "")
      write.phy(a[[i]], paste(gene, b, "masked.phy", sep = "-"))
      rownames(a[[i]]) <- gsub("-", "_", rownames(a[[i]]))
      write.nex(a[[i]], paste(gene, b, "masked.nex", sep = "-"))
    }
    a <- do.call(cbind.DNAbin, c(a, fill.with.gaps = TRUE))
    b <- paste(length(n), "blocks", sep = "")
    write.phy(a, paste(gene, "masked.phy", sep = "-"))
    rownames(a) <- gsub("-", "_", rownames(a))
    write.nex(a, paste(gene, "masked.nex", sep = "-"))
  }
  
  ## summary of result
  ## -----------------
#   slog(paste("\n\n--- final alignment of", gene, "---"),
#        paste("\nnumber of sequences     :", nrow(a)),
#        paste("\nnumber of base pairs    :", ncol(a)),
#        file = logfile)
#   if ( length(n) > 1 ) {
#     slog("\nnumber of blocks        :", length(n), 
#          "(", paste(n, collapse = ", "), ")", 
#          file = logfile) 
#   }
  
  dbDisconnect(conn)
  slog("\n\nSTEP I finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
  invisible(x)
}
