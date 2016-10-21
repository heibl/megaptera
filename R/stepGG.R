## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-08-16)

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
  logfile <- paste(gene, "stepG.log", sep = "-")
  if ( !quiet & file.exists(logfile) ) unlink(logfile)
  if ( !quiet )  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
                      paste("\n", Sys.time(), sep = ""),
                      "\nSTEP G: alignment\n", 
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
  tax <- dbReadTaxonomy(x, subset = msa.tab)
  gt <- tax2tree(tax, tip.rank = "spec")
  gt <- decomposePhylo(gt)
  subtree <- sort(unique(gt$subtree.set$subtree))
  
  ## alignment of genera/subtrees -- either sequential or parallel
  ## ----------------------------------------------------
  SQL <- paste("SELECT status = 'subtree-aligned' AS status FROM", msa.tab)
  if ( !all(dbGetQuery(conn, SQL)$status) ){
    
    slog("\n.. aligning", length(subtree),
         "subtrees ..", file = logfile)
    cpus <- x@params@cpus
    if ( length(gt$subtree.set) < cpus | !x@params@parallel ){
      lapply(subtree, alignSubtree, 
             megProj = x, taxon = gt$subtree.set)
      cluster.open <- FALSE
    } else {
      slog("\n", file = logfile)
      sfInit(parallel = TRUE, cpus = cpus, 
             type = x@params@cluster.type)
      sfLibrary("megaptera", character.only = TRUE)
      sfExport("gt$subtree.set", "x")
      seqs <- sfLapply(gt$subtree.set, alignSubtree, megProj = x)
      cluster.open <- TRUE
    }
  } else {
    cluster.open <- FALSE
  }
  ## if dataset contains only one genus, 
  ## there is no need
  ## for alignment along guide tree
  ## ------------------------------
  # if ( length(seqs) == 1 ){
  #   seqs <- seqs[[1]]
  # } else {
  
  ## prepare guide tree for alignment
  ## --------------------------------
  # if ( inherits(x@taxon, "taxonGuidetree") ){
  #   gt <- comprehensiveGuidetree(x, tip.rank = "gen")
  # } else { 
  #   gt <- tax2tree(tax, tip.rank = "gen")
  # }
  # gt <- drop.tip(gt, pruned <- setdiff(gt$tip.label, 
  #                                      names(seqs)))
  # if ( !is.binary.tree(gt) ){
  #   gt <- multi2di(gt)
  # }
  # 
  # if ( length(pruned) > 0 ) 
  #   slog("\n.. pruning", length(pruned), 
  #        "taxa from guide tree ..", file = logfile)
  # seqs <- seqs[intersect(names(seqs), gt$tip.label)]
  
  ## loop over internal nodes of guide tree
  ## --------------------------------------
  gt <- gt$guidetree
  slog("\n.. aligning sister clades higher than genera ..", 
       file = logfile)
  i <- 1
  # while ( Ntip(gt) >= 2 ){
          for ( i in 2:13 ){
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
        sfExport("tp", "seqs") # megProj already exported
        xx <- sfSapply(tp, alignSisterClades, seqs = seqs, 
                       megProj = x)
      } else {
        slog(" serially", file = logfile)
        xx <- sapply(tp, mergeSubMSA, x = x)
      }
    }
    
    # if ( Nnode(gt) == 1 ) { ## basale Polytomie
    #   seqs <- xx
    #   i <- i + 1
    #   break
    # }
    keep.tips <- sapply(tp, head, 1)
    drop.tips <- unlist(lapply(tp, tail, -1))
    gt$tip.label[match(keep.tips, gt$tip.label)] <- xx
    gt <- drop.tip(gt, drop.tips)
    i <- i + 1
  } # end of WHILE-loop
  
  ## align outgroup
  # if ( !quiet ) slog("\nLEVEL", i, file = logfile)
  # if ( length(seqs) > 2 )
  #   stop("uncomplete alignment in WHILE loop")
  # if ( length(seqs) == 2 )
  # seqs <- mergeSubMSA(x, )
  # else seqs <- seqs[[1]]
  # }
  
  ## prune ends of alignment from 'thin tails'
  ## -----------------------------------------
  seqs <- trimEnds(seqs, min(nrow(seqs), min.n.seq))
  seqs <- deleteEmptyCells(seqs, quiet = TRUE)
  
  ## sort alignment taxonomically
  ## ----------------------------
  if ( nrow(seqs) > 1 ){
    rownames(seqs) <- gsub("_R_", "", rownames(seqs))
    gt <- tax2tree(tax, tip.rank = tip.rank)
    gt <- ladderize(gt)
    gt <- drop.tip(gt, setdiff(gt$tip.label, rownames(seqs)))
    seqs <- seqs[match(gt$tip.label, rownames(seqs)), ]
  }
  
  ## write to database -- serially or in parallel
  ## --------------------------------------------
  slog("\n.. write alignment to database ..", file = logfile)
  if ( cluster.open ){
    if ( nrow(seqs) > 1000 ){
      id <- seq(from = 1, to = nrow(seqs), by = ceiling(nrow(seqs)/x@params@cpus))
      id <- data.frame(from = id, to = c(id[-1] - 1, nrow(seqs)))
      aa <- apply(id, 1, function(i, a) a[i[1]:i[2], ], a = seqs)
      sfLibrary("seqinr", character.only = TRUE)
      sfExport("aa")
      sfLapply(x = aa, fun = dbWriteMSA, megapteraProj = x, 
               status = "aligned")
    } else {
      dbWriteMSA(x, dna = seqs, status = "aligned")
    }
    sfStop()
  } else {
    dbWriteMSA(x, dna = seqs, status = "aligned")
  }
  ## This is just a hack to allow for negative matching of regexs
  ## In the future stepG should insert the highest non-saturated taxonomic rank
  dbSendQuery(conn, paste("UPDATE", msa.tab, 
                          "SET status = 'included'", 
                          "WHERE status is NULL"))
  dbDisconnect(conn)
  
  ## write files
  ## -----------
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
      "\n\nSTEP G finished",
      file = logfile)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), file = logfile)
  }
  invisible(x)
}
