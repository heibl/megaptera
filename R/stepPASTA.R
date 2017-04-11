## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-03-28)

#' @export
#' @import DBI snow snowfall
#' @importFrom snow setDefaultClusterOptions

stepPASTA <- function(x, k = 200){	
  
  start <- Sys.time()
  quiet = FALSE
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") 
    stop("undefined locus not allowed")
  
  ## check if previous step has been run
  ## -----------------------------------
  status <- dbProgress(x)
  if (status$step_f == "pending") {
    stop("the previous step has not been called yet")
  }
  if (status$step_f == "error") {
    stop("the previous step has terminated with an error")
  }
  if (status$step_f == "failure") {
    slog("\nNo data from upstream available - quitting", file = "")
    dbProgress(x, "step_g", "failure")
    return()
  }
  if (status$step_f == "success") {
    dbProgress(x, "step_g", "error")
  }
  
  ## PARAMETERS
  ## ----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- match.arg(x@taxon@tip.rank, c("species", "genus"))
  msa.tab <- paste(x@taxon@tip.rank, gsub("^_", "", gene), sep = "_")
  min.n.seq <- x@params@min.n.seq
  align.exe <- x@align.exe
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepPASTA.log")
  if (!quiet & file.exists(logfile)) unlink(logfile)
  if (!quiet)  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
                      paste("\n", Sys.time(), sep = ""),
                      "\nSTEP G: PASTA alignment\n", 
                      paste("\n.. locus:", x@locus@sql), file = logfile)
  
  ## open database connection
  conn <- dbconnect(x)
  
  ## check if msa table exists
  ## -------------------------
  if (!dbExistsTable(conn, msa.tab)){
    dbDisconnect(conn)
    slog("\nWARNING: table", msa.tab, "does not exist!\n", 
         file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP G finished after", round(td, 2), attr(td, "units"), 
         "\n", file = logfile)
    return()
  }
  
  ## check if at least 3 (ingroup) species are available
  ## ---------------------------------------------------
  n <- paste("SELECT count(taxon) FROM", msa.tab)
  n <- dbGetQuery(conn, n)$count
  if (n < 3){
    dbDisconnect(conn)
    slog("\nWARNING: only", n, "species available, no alignment possible\n", 
         file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP GG finished after", round(td, 2), attr(td, "units"), 
         "\n", file = logfile)
    dbProgress(x, "step_g", "failure")
    return()
  }
  if (n < 100){ # 100 is arbitrary
    n <- paste("SELECT taxon FROM", msa.tab)
    n <- dbGetQuery(conn, n)
    n <- which(is.ingroup(x, n$taxon))
    if (length(n) < 3){
      dbDisconnect(conn)
      slog("\nWARNING:", length(n), "ingroup species available, no alignment possible\n", 
           file = logfile)
      td <- Sys.time() - start
      slog("\nSTEP GG finished after", round(td, 2), attr(td, "units"), 
           "\n", file = logfile)
      dbProgress(x, "step_g", "failure")
      return()
    }
  }
  
  ## read taxonomy and create guide tree
  ## -----------------------------------
  gt <- comprehensiveGuidetree(x, tip.rank = tip.rank, subset = msa.tab)
  guidetree <- gt
  
  ## less then k species will be aligned with MAFFT-LINSI
  ## ----------------------------------------------------
  if (Ntip(gt) <= k){
    slog("\n..", Ntip(gt), "species will be aligned with MAFFT L-INS-i", 
         file = logfile)
    seqs <- dbReadDNA(x, msa.tab)
    if (is.matrix(seqs)) seqs <- del.gaps(seqs)
    seqs <- mafft(seqs, method = "localpair", gt = gt, exec = align.exe)
    
  ## more than k species will be aligned with PASTA
  ## ----------------------------------------------
  } else {
    
    gt <- decomposePhylo(gt)
    subtree <- sort(unique(gt$subtree.set$subtree))
    
    ## alignment of subsets -- either sequential or parallel
    ## -----------------------------------------------------
    SQL <- paste("SELECT status = 'subtree-aligned' AS status FROM", msa.tab)
    if (!all(dbGetQuery(conn, SQL)$status)){
      
      slog("\n.. aligning", length(subtree),
           "subtrees ..", file = logfile)
      cpus <- x@params@cpus
      if (Ntip(gt$guidetree) < cpus | !x@params@parallel){
        lapply(subtree, alignSubtree, 
               megProj = x, taxon = gt$subtree.set)
        cluster.open <- FALSE
      } else {
        slog("\n", file = logfile)
        sfInit(parallel = TRUE, cpus = cpus, 
               type = x@params@cluster.type)
        sfLibrary("megaptera", character.only = TRUE)
        sfExport("gt", "x")
        sfLapply(subtree, alignSubtree, megProj = x, taxon = gt$subtree.set)
        sfStop()
        cluster.open <- TRUE
      }
    } else {
      cluster.open <- FALSE
    }
    
    ## Step 2: calculate minimum spanning tree
    ## ---------------------------------------
    gt <- phylo2mst(gt$guidetree)
    
    ## Step 3: merge "Type 1 sub-alignments"
    ## -------------------------------------
    
    ## get names of neighboring pairs:
    gt[upper.tri(gt)] <- 0
    id <- which(gt == 1, arr.ind = TRUE)
    np <- data.frame(a = rownames(gt)[id[, 1]],
                     b = colnames(gt)[id[, 2]])
    np <- np[order(np$a), ]
    np <- apply(np, 1, as.list)
    np <- lapply(np, unlist)
    
    ## merge pairs of MSAs
    if (!cluster.open){
      slog(" serially", file = logfile)
      seqs <- lapply(np, mergeSubMSA2, megProj = x)
    } else {
      if (length(np) > 1){
        slog(" in parallel on", min(nrow(np), x@params@cpus), 
             "nodes", file = logfile)
        sfInit(parallel = TRUE, cpus = x@params@cpus, 
               type = x@params@cluster.type)
        sfLibrary("megaptera", character.only = TRUE)
        sfExport(list = c("x", "np"), debug = TRUE)
        # sfExport("np") # x already exported
        seqs <- sfLapply(np, mergeSubMSA2, megProj = x)
        sfStop()
      } else {
        slog(" serially", file = logfile)
        seqs <- sapply(np, mergeSubMSA2, x = x)
      }
    }
    
    ## Step 4: TRANSIVITY MERGE
    ## ------------------------
    
    ## A: parallel:
    
    sfInit(parallel = TRUE, cpus = x@params@cpus,
           type = x@params@cluster.type)
    sfLibrary("megaptera", character.only = TRUE)
    
    ## alignments as list of indices
    seqs <- lapply(seqs, DNAbin2index)
    
    repeat{
      ## create list of overlapping pairs
      op <- overlappingPairs2(seqs)
      
      ## transitivity merger
      sfExport("op", "seqs")
      seqs2 <- sfLapply(op, transitivityMerge, obj = seqs)
      seqs <- seqs2
      # seqs2 <- lapply(op, transitivityMerge, obj = seqs)
      if ( length(seqs) == 1 ) break
    }
    
    sfStop()
    
    ## B: serially
    ## -----------
    ## alignments as list of indices
    # seqs <- lapply(seqs, DNAbin2index)
    # 
    # ## very rough loop!
    # repeat {
    #   s <- names(seqs[[1]])
    #   s <- sapply(seqs, 
    #               function(z, s) length(intersect(names(z), s)) > 0, 
    #               s = s)
    #   s <- which(s)[1:2]
    #   
    #   seqs <- c(seqs, list(transitivityMerge(obj = seqs, s)))
    #   seqs <- seqs[-s]
    #   cat("\nnumber of subMSAs:", length(seqs))
    #   if ( length(seqs) == 1 ) break
    # }
    seqs <- seqs[[1]]
    
    dna <- dbReadDNA(x, msa.tab) 
    dna <- del.miss(dna)
    
    seqs <- index2DNAbin(dna, seqs)
  }

  
  ## prune ends of alignment from 'thin tails'
  ## -----------------------------------------
  seqs <- trimEnds(seqs, min(nrow(seqs), min.n.seq))
  seqs <- deleteEmptyCells(seqs, quiet = TRUE)

  ## sort alignment taxonomically
  ## ----------------------------
  if (nrow(seqs)){
    # rownames(seqs) <- gsub("_R_", "", rownames(seqs))
    gt <- ladderize(guidetree)
    gt <- drop.tip(gt, setdiff(gt$tip.label, rownames(seqs)))
    seqs <- seqs[match(gt$tip.label, rownames(seqs)), ]
  }
  
  ## write to database 
  ## -----------------
  dbWriteMSA(x, dna = seqs, status = "aligned")
  dbDisconnect(conn)
  
  # write files
  # -----------
  write.phy(seqs, paste0("msa/", gene, ".phy"))
  rownames(seqs) <- gsub("-", "_", rownames(seqs))
  write.nex(seqs, paste0("msa/", gene, ".nex"))

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
    slog(" after", round(td, 2), attr(td, "units"), "\n", file = logfile)
  }
  dbProgress(x, "step_g", "success")
}
