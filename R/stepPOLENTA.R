## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2018-12-18)

#' @title STEP G: POLENTA Alignment
#' @description Use POLENTA to align tip-rank-level sequences.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param k An integer giving the size of a subset.
#' @note This function is not fully implemented!
#' @import DBI snow snowfall
#' @importFrom snow setDefaultClusterOptions
#' @export

stepPOLENTA <- function(x, k = 200){	
  
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
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepPASTA.log")
  if (!quiet & file.exists(logfile)) unlink(logfile)
  if (!quiet) slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
                   paste("\n", Sys.time()),
                   "\nSTEP G: POLENTA alignment\n", 
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
  n <- paste("SELECT count(taxon)",
             "FROM", msa.tab,
             "WHERE", wrapSQL(gene, "locus", "="))
  n <- dbGetQuery(conn, n)$count
  if (n < 3){
    dbDisconnect(conn)
    slog("\nWARNING: only", n, "species available, no alignment possible\n", 
         file = logfile, megProj = x)
    td <- Sys.time() - start
    slog("\nSTEP GG finished after", round(td, 2), attr(td, "units"), 
         "\n", file = logfile, megProj = x)
    dbProgress(x, "step_g", "failure")
    return()
  }
  if (n < 100){ # 100 is arbitrary
    n <- dbGetQuery(conn, paste("SELECT taxon FROM", msa.tab,
                                "WHERE", wrapSQL(gene, "locus", "=")))
    n <- which(is.ingroup(x, n$taxon))
    if (length(n) < 3){
      dbDisconnect(conn)
      slog("\nWARNING:", length(n), "ingroup species available, no alignment possible\n", 
           file = logfile, megProj = x)
      td <- Sys.time() - start
      slog("\nSTEP GG finished after", round(td, 2), attr(td, "units"), 
           "\n", file = logfile, megProj = x)
      dbProgress(x, "step_g", "failure")
      return()
    }
  }
  
  ## read taxonomy and create guide tree
  ## -----------------------------------
  gt <- comprehensiveGuidetree(x, tip.rank = tip.rank, subset = msa.tab)
  guidetree <- gt
  
  ## POLENTA
  # seqs <- polenta(seqs, k = k, parallel = x@params@parallel,
  #                 exec = x@align.exe, ncore = x@params@cpus)

  
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
  if (!quiet) {
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
