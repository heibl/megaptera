## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-02-20)

#' @export
#' @import DBI
#' @import snowfall
#' @importFrom snowfall sfInit

stepGG <- function(x){	
  
  start <- Sys.time()
  quiet = FALSE
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  if ( x@locus@kind == "undefined" ) 
    stop("undefined locus not allowed")
  STATUS <- checkStatus(x)
  if ( !all(STATUS[1:6]) ){
    stop("step", names(STATUS)[min(which(!STATUS))] ,
         " has not been run")
  }
  
  ## PARAMETERS
  ## ----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- x@taxon@tip.rank
  msa.tab <- paste(tip.rank, gsub("^_", "", gene), sep = "_")
  min.n.seq <- x@params@min.n.seq
  align.exe <- x@align.exe
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepGG.log")
  if ( !quiet & file.exists(logfile) ) unlink(logfile)
  if ( !quiet )  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
                      paste("\n", Sys.time(), sep = ""),
                      "\nSTEP GG: alignment\n", 
                      paste("\n.. locus:", x@locus@sql), file = logfile)
  
  
  
  ## open database connection
  conn <- dbconnect(x)
  
  ## check if msa table exists
  ## -------------------------
  if ( !dbExistsTable(conn, msa.tab) ){
    dbDisconnect(conn)
    slog("\nWARNING: table", msa.tab, "does not exist!\n", 
         file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP G finished after", round(td, 2), attr(td, "units"), 
         "\n", file = logfile)
    return()
  }
  
  ## check if at least 3 species are available
  ## -----------------------------------------
  n <- paste("SELECT count(spec) FROM", msa.tab)
  n <- dbGetQuery(conn, n)$count
  if ( n < 3 ){
    dbDisconnect(conn)
    slog("\nWARNING: only", n, "species available, no alignment possible\n", 
         file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP GG finished after", round(td, 2), attr(td, "units"), 
         "\n", file = logfile)
    return()
  }
  
  ## Check if alignment is necessary.
  ## Any 'raw' entry in the status column
  ## is a trigger for aligning/updating; 
  ## 'raw' is inserted by stepF and dbExcludeSpec.
  ## --------------------------------------------
  # SQL <- paste("SELECT status = 'raw' AS status FROM", msa.tab)
  # if ( !any(dbGetQuery(conn, SQL)$status) ){
  #   dbDisconnect(conn)
  #   slog(paste("\ntable '", msa.tab, "' is up to date", 
  #              " - nothing to be done\n", sep = ""), 
  #        file = logfile)
  #   td <- Sys.time() - start
  #   slog("\nSTEP G finished after", round(td, 2), attr(td, "units"), 
  #        "\n", file = logfile)
  #   return()
  # }
  
  ## read taxonomy relation from database
  ## ------------------------------------
  if ( inherits(x@taxon, "taxonGuidetree") ){
    gt <- comprehensiveGuidetree(x, tip.rank = "spec", subset = msa.tab)
  } else { 
    tax <- dbReadTaxonomy(x, subset = msa.tab)
    gt <- tax2tree(tax, tip.rank = "spec")
  }
  guidetree <- gt
  gt <- decomposePhylo(gt)
  subtree <- sort(unique(gt$subtree.set$subtree))
  
  ## alignment of genera/subtrees -- either sequential or parallel
  ## ----------------------------------------------------
  SQL <- paste("SELECT status = 'subtree-aligned' AS status FROM", msa.tab)
  if ( !all(dbGetQuery(conn, SQL)$status) ){
    
    slog("\n.. aligning", length(subtree),
         "subtrees ..", file = logfile)
    cpus <- x@params@cpus
    if ( Ntip(gt$guidetree) < cpus | !x@params@parallel ){
      lapply(subtree, alignSubtree, 
             megProj = x, taxon = gt$subtree.set)
      cluster.open <- FALSE
    } else {
      slog("\n", file = logfile)
      sfInit(parallel = TRUE, cpus = cpus, 
             type = x@params@cluster.type)
      sfLibrary("megaptera", character.only = TRUE)
      sfExport("gt", "x")
      seqs <- sfLapply(subtree, alignSubtree, 
                       megProj = x, taxon = gt$subtree.set)
      cluster.open <- TRUE
    }
  } else {
    cluster.open <- FALSE
  }
  
  ## loop over internal nodes of guide tree
  ## --------------------------------------
  gt <- gt$guidetree
  slog("\n.. aligning sister clades higher than genera ..", 
       file = logfile)
  i <- 1
  while ( Ntip(gt) >= 2 ){
          # for ( i in 2:48  ){
    #         load("~/r/dicaryota/data/BUGSEARCH.rda.RData")
    if ( !quiet ) slog("\n\nLEVEL", i, file = logfile)
    ## tp: list of 'terminal clades', i.e. clades that
    ## consist only of tips
    tp <- terminalSisters(gt)
    ntp <- length(tp)
    slog(": aligning", ntp, ifelse(ntp == 1, "pair", "pairs"),
         "of sister clades", file = logfile)
    
    ## align each set of sequences in list 'seqs'
    ## -- either serially or in parallel
    ## ---------------------------------
    if ( !cluster.open ){
      slog(" serially", file = logfile)
      xx <- sapply(tp, mergeSubMSA, x = x)
    } else {
      if ( ntp > 1 ){
        slog(" in parallel on", min(ntp, x@params@cpus), 
             "nodes", file = logfile)
        sfLibrary("megaptera", character.only = TRUE)
        sfExport("tp") # x already exported
        xx <- sfSapply(tp, mergeSubMSA, x = x)
      } else {
        slog(" serially", file = logfile)
        xx <- sapply(tp, mergeSubMSA, x = x)
      }
    }
    keep.tips <- sapply(tp, head, 1)
    drop.tips <- unlist(lapply(tp, tail, -1))
    gt$tip.label[match(keep.tips, gt$tip.label)] <- xx
    if ( Ntip(gt) == 2 ) break
    gt <- drop.tip(gt, drop.tips)
    i <- i + 1
  } # end of WHILE-loop
  
  seqs <- dbReadDNA(x, msa.tab, regex = TRUE, blocks = "ignore")
  
  ## prune ends of alignment from 'thin tails'
  ## -----------------------------------------
  seqs <- trimEnds(seqs, min(nrow(seqs), min.n.seq))
  seqs <- deleteEmptyCells(seqs, quiet = TRUE)

  ## sort alignment taxonomically
  ## ----------------------------
  if ( nrow(seqs) > 1 ){
    # rownames(seqs) <- gsub("_R_", "", rownames(seqs))
    gt <- ladderize(guidetree)
    gt <- drop.tip(gt, setdiff(gt$tip.label, rownames(seqs)))
    seqs <- seqs[match(gt$tip.label, rownames(seqs)), ]
  }

  dbDisconnect(conn)
  
  # write files
  # -----------
  write.phy(seqs, paste(gene, "phy", sep = "."))
  rownames(seqs) <- gsub("-", "_", rownames(seqs))
  write.nex(seqs, paste(gene, "nex", sep = "."))

  ## calculate mean absolute deviation (MAD)
  ## see Smith, Beaulieau, Donoghue (2009)
  this.mad <- round(MAD(seqs), 5)
  slog("\n.. mean absolute deviation:",
       this.mad, "..", file = logfile)

  ## summary
  ## -------
  if ( !quiet ) {
    slog(
      paste("\n\n--- final alignment of", gene, "---"),
      paste("\nnumber of sequences     :", nrow(seqs)),
      paste("\nnumber of base pairs    :", ncol(seqs)),
      paste("\nmean absolute deviation :", this.mad),
      "\n\nSTEP GG finished",
      file = logfile)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), file = logfile)
  }
  invisible(x)
}
