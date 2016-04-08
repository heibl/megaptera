## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-04-08)

stepH <- function(x, max.mad = 0.01){
  
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
  min.n.seq <- x@params@min.n.seq
  max.bp <- x@params@max.bp * 1.5
  
  ## iniate logfile
  ## --------------
  logfile <- paste(gene, "stepH.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP H: detection and separation of unsaturated blocks\n", 
       paste("\n.. locus:", x@locus@sql), file = logfile)
  
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
    slog("\nSTEP H finished after", round(td, 2), attr(td, "units"), "\n",
         file = logfile)
    return()
  }
  
  ## read alignment
  ## --------------
  slog("\n.. reading alignment with ", file = logfile)
  if ( x@params@parallel ){
    spec <- dbGetQuery(conn, paste("SELECT", tip.rank, "FROM", msa.tab))[, tip.rank]
    spec <- paste("^", spec, "$", sep = "")
    n <- length(spec)
    id <- seq(from = 1, to = n, by = ceiling(n/x@params@cpus))
    id <- data.frame(from = id, to = c(id[-1] - 1, n))
    spec <- apply(id, 1, function(i, a) paste(a[i[1]:i[2]], collapse = "|"), a = spec)
    sfInit(parallel = TRUE, cpus = x@params@cpus, 
           type = x@params@cluster.type)
    sfLibrary("megaptera", character.only = TRUE)
    sfLibrary("seqinr", character.only = TRUE)
    megProj <- x
    sfExport("spec", "msa.tab", "megProj")
    a <- sfLapply(x = megProj, fun = dbPReadDNA, tab.name = msa.tab, 
                  regex = TRUE, max.bp = max.bp)
    sfStop()
    a <- do.call(rbind, a)
  } else {
    a <- dbReadDNA(x, msa.tab, taxon = ".+", regex = TRUE,
                   ignore.excluded = TRUE, blocks = "ignore",
                   max.bp = max.bp)
  }
  ## check if stepG has been run properly
  ## should include checking for status = 'raw'
  if ( !is.matrix(a) ){
    stop("stepG has not been called yet")
  }
  ## how many species are in the table/alignment
  slog(nrow(a), "species", file = logfile)
  
  ## special case: alignment consists of one species
  ## SHOULD BE MADE IMPOSSIBLE IN stepG
  if ( nrow(a) == 1 ) {
    dbDisconnect(conn)
    slog("\n", file = logfile)
    dbWriteMSA(x, dna = a, masked = FALSE)
    return()
  }
  
  ## Set status column to 'aligned':
  ## If this is not done and the alignment is below
  ## MAD threshold, blocks from previous runs will
  ## not be eliminated.
  ## -----------------
  SQL <- paste("UPDATE", msa.tab,
               "SET status = 'aligned'",
               "WHERE status ~ 'block'")
  dbSendQuery(conn, SQL)
  
  ## checking saturation of alignment as measures
  ## by the median average difference (MAD) between
  ## corrected and JC69-corrected patristic distances
  ## ------------------------------------------------
  check.mad <- MAD(a)
  slog("\n.. assessing saturation: MAD =", round(check.mad, 5),
       file = logfile)
  if ( check.mad <= max.mad ){
    dbDisconnect(conn)
    slog(paste("\n.. MAD below threshold of ", max.mad, ": ",
         "alignment not saturated\n", sep = ""), file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP H finished after", round(td, 2), attr(td, "units"), "\n",
         file = logfile)
    return()
  }
  
  ## prepare guide tree for alignment
  ## --------------------------------
  slog("\n.. preparing guidetree", file = logfile)
  if ( inherits(x@taxon, "taxonGuidetree") ){
    gt <- comprehensiveGuidetree(x, tip.rank = tip.rank, subset = a)
  } else { 
    gt <- dbReadTaxonomy(x, subset = a)
    gt <- tax2tree(gt, tip.rank = tip.rank)
  }
  gt <- fixNodes(gt)
  a <- a[match(gt$tip.label, rownames(a)), ]
  a <- list(a)
  this.root <- Ntip(gt) + 1
  
  slog("\n.. successively splitting up alignment", file = logfile)
  repeat {

    aa <- splitAlignment(this.root, gt, a[[1]])
    a <- c(a, aa)
    a <- a[-1]
    check.mad <- sapply(a, MAD)
    slog("\n   - maximum MAD:", max(check.mad), file = logfile)
    if ( all(check.mad <= max.mad) ) break
    a <- a[order(check.mad, decreasing = TRUE)]
    this.root <- noi(gt, rownames(a[[1]]))
    
  }
  slog("\n.. alignment was split into", length(a), 
       "unsaturated blocks", file = logfile)
  check.size <- sapply(a, nrow)
  a <- a[order(check.size, decreasing = TRUE)]
  
  ## write alignment blocks to database
  ## ----------------------------------
  for ( i in 1:length(a) ){
    #     dbWriteMSA(megapteraProj = x, 
    #                dna = a[[i]], 
    #                status = paste("block", i), 
    #                masked = FALSE)
    tips <- lapply(a, rownames)
    names(tips) <- paste("block", 1:length(tips))
    tips <- sapply(tips, wrapSQL, term = tip.rank, boolean = "OR", operator = "=")
    SQL <- paste("UPDATE", msa.tab,
                 "SET", wrapSQL(names(tips), term = "status", 
                                boolean = NULL, operator = "="),
                 "WHERE", tips)
    lapply(SQL, dbSendQuery, conn = conn)
  }
  
  ## write files
  ## -----------
  #   slog("\n.. write alignment to files ..", file = logfile)
  #   write.phy(a, paste(gene, "masked.phy", sep = "-"))
  #   rownames(a) <- gsub("-", "_", rownames(a))
  #   write.nex(a, paste(gene, "masked.nex", sep = "-"))
  
  ## summary of result
  ## -----------------
  #     slog(paste("\n\n--- final alignment of", gene, "---"),
  #          paste("\nnumber of sequences     :", nrow(a)),
  #          paste("\nnumber of base pairs    :", ncol(a)),
  #          file = logfile)
  #     if ( length(blk) > 1 ) {
  #       slog("\nnumber of blocks        :", length(blk), 
  #            "(", paste(blk, collapse = ", "), ")", 
  #            file = logfile) 
  #     }
  #     slog(paste("\nmean absolute deviation :", this.mad), 
  #          file = logfile)
  #   }
  dbDisconnect(conn)
  slog("\n\nSTEP H finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
}
