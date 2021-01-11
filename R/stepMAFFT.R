## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2019-10-30)

#' @title STEP G: MAFFT Alignment
#' @description Use MAFFT to align tip-rank-level sequences.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param method A character string giving the alignment method. Available
#'   accuracy-oriented methods for less than 200 sequences are
#'   \code{"localpair"}, \code{"globalpair"}, and \code{"genafpair"};
#'   \code{"retree 1"} and \code{"retree 2"} are for speed-oriented alignment.
#'   The default is \code{"auto"}, which lets MAFFT choose an appropriate
#'   alignment method.
#' @param maxiterate An integer giving the number of cycles of iterative
#'   refinement to perform. Possible choices are \code{0}: progressive method,
#'   no iterative refinement (default); \code{2}: two cycles of iterative
#'   refinement; \code{1000}: at most 1000 cycles of iterative refinement.
#' @param op A numeric giving the \code{gap opening penalty} at group-to-group
#'   alignment; default 1.53.
#' @param ep A numeric giving the offset value, which works like \code{gap
#'   extension penalty}, for group-to-group alignment; default 0.0, but 0.123 is
#'   recommended if no long indels are expected.
#' @details \code{"localpair"} selects the \bold{L-INS-i} algorithm, probably
#'   most accurate; recommended for <200 sequences; iterative refinement method
#'   incorporating local pairwise alignment information.
#'
#'   \code{"globalpair"} selects the \bold{G-INS-i} algorithm suitable for
#'   sequences of similar lengths; recommended for <200 sequences; iterative
#'   refinement method incorporating global pairwise alignment information.
#'
#'   \code{"genafpair"} selects the \bold{E-INS-i} algorithm suitable for
#'   sequences containing large unalignable regions; recommended for <200
#'   sequences.
#'
#'   \code{"retree 1"} selects the \bold{FFT-NS-1} algorithm, the simplest
#'   progressive option in MAFFT; recommended for >200 sequences.
#'
#'   \code{"retree 2"} selects the \bold{FFT-NS-2} algorithm that uses a second
#'   iteration of alignment based on a guide tree computed from an FFT-NS-1
#'   aligment; this is the default in MAFFT; recommended for >200 sequences.
#' @return None, \code{stepMAFFT} is called for its side effect.
#' @references Katoh, K. and H. Toh. 2008. Recent developments in the MAFFT
#'   multiple sequence alignment program. \emph{Briefings in Bioinformatics}
#'   \bold{9}: 286-298.
#'
#'   Katoh, K., K.-i. Kuma, H. Toh, and T. Miyata. 2005. Mafft version 5:
#'   improvement in accuracy of multiple sequence alignment. \emph{Nucleic Acids
#'   Research} \bold{33}: 511--518.
#'
#'   Katoh, K., K. Misawa, K.-i. Kuma, and T. Miyata. 2002. Mafft: a novel
#'   method for rapid multiple sequence alignment based on fast Fourier
#'   transform. \emph{Nucleid Acids Research} \bold{30}: 3059--3066.
#'
#'   \url{http://mafft.cbrc.jp/alignment/software/index.html}
#' @importFrom ape ladderize
#' @importFrom DBI dbDisconnect dbExistsTable
#' @importFrom ips write.nex
#' @export

stepMAFFT <- function(x, method = "auto", maxiterate = 0, op = 1.53, ep = 0){	
  
  start <- Sys.time()
  quiet <- FALSE
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") 
    stop("undefined locus not allowed")
  
  ## check if previous step has been run
  ## -----------------------------------
  # status <- dbProgress(x)
  # if (status$step_f == "pending") {
  #   stop("the previous step has not been called yet")
  # }
  # if (status$step_f == "error") {
  #   stop("the previous step has terminated with an error")
  # }
  # if (status$step_f == "failure") {
  #   slog("\nNo data from upstream available - quitting", file = "")
  #   dbProgress(x, "step_g", "failure")
  #   return()
  # }
  # if (status$step_f == "success") {
  #   dbProgress(x, "step_g", "error")
  # }
  
  ## PARAMETERS
  ## ----------
  gene <- x@locus@sql
  tip.rank <- match.arg(x@taxon@tip.rank, c("species", "genus"))
  min.n.seq <- x@params@min.n.seq
  dbl <- x@params@debug.level
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepPASTA.log")
  if (!quiet & file.exists(logfile)) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
                    paste0("\n", Sys.time()),
                    "\nSTEP G: MAFFT alignment\n", 
                    paste("\nLocus:", gene), file = logfile, megProj = x)
  
  ## open database connection
  conn <- dbconnect(x)
  
  ## check if at least 3 (ingroup) species are available
  ## ---------------------------------------------------
  n <- paste("SELECT count(taxon)",
             "FROM sequence_selected",
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
    n <- dbGetQuery(conn, paste("SELECT taxon FROM sequence_selected", 
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
  
  ## Read taxonomy and create guide tree
  ## -----------------------------------
  slog("\nCreating comprehensive guide tree ... ", file = logfile, megProj = x)
  gt <- comprehensiveGuidetree(x, tip.rank = tip.rank, subset = "sequence_selected")
  slog("OK", file = logfile, megProj = x)
  
  ## Read DNA sequences and check if guide tree is compatible
  ## --------------------------------------------------------
  slog("\nReading sequences ... ", file = logfile, megProj = x)
  seqs <- dbReadMSA(x, label = "taxon")
  slog("OK", file = logfile, megProj = x)
  
  ## Check ...
  if (is.matrix(seqs)){
    seq_names <- rownames(seqs)
  } else {
    seq_names <- names(seqs)
  }
  ## 1. if alignment
  not_in_seqs <- setdiff(gt$tip.label, seq_names)
  not_in_gt <- setdiff(seq_names, gt$tip.label)
  if (length(not_in_gt) | length(not_in_seqs)){
    if (dbl > 0){
      if (length(not_in_gt)){
        slog("\n", length(not_in_gt), " sequence(s) missing in guide tree:\n- ", 
             paste(head(not_in_gt), collapse = "\n- "), 
             file = logfile, sep = "", megProj = x)
      }
      if (length(not_in_seqs)){
        slog("\n", length(not_in_seqs), " tip(s) missing from sequences:\n- ", 
             paste(head(not_in_seqs), collapse = "\n- "), 
             file = logfile, sep = "", megProj = x)
      }
      if (dbl > 3){
        err_file <- paste0("log/", Sys.Date(), "_DEBUG_stepMafft.rda")
        save(x, seqs, gt, file = err_file)
        slog("\nSequences and guide tree written to '", err_file, "'\n", 
             sep = "", file = logfile)
      }
    }
    stop("sequences names do not match names in guide tree")
  }
  
  ## MAFFT alignment
  ## ---------------
  slog("\nAligning sequences with MAFFT ... ", file = logfile, megProj = x)
  seqs <- mafft(seqs, method = method, gt = gt, 
                maxiterate = maxiterate, op = op, ep = ep,
                thread = x@params@cpus, exec = x@align.exe)
  if (is.logical(seqs)){
    slog("failed", file = logfile, megProj = x)
    stop("MAFFT alignment failed")
  }
  slog("OK", file = logfile, megProj = x)
  
  ## Prune ends of alignment from 'thin tails'
  ## -----------------------------------------
  slog("\nPruning 'thin tails' from alignment ... ", file = logfile, megProj = x)
  seqs <- trimEnds(seqs, min(nrow(seqs), min.n.seq))
  slog("OK", file = logfile, megProj = x)
  
  ## Detect and delete empty cells from alignment
  ## --------------------------------------------
  slog("\nDetect and delete empty cells from alignment ... ", file = logfile, megProj = x)
  seqs <- deleteEmptyCells(seqs, quiet = TRUE)
  slog("OK", file = logfile, megProj = x)

  ## Sort alignment taxonomically
  ## ----------------------------
  if (nrow(seqs)){
    slog("\nSort alignment according to guide tree ... ", file = logfile, megProj = x)
    gt <- ladderize(gt)
    gt <- drop.tip(gt, setdiff(gt$tip.label, rownames(seqs)))
    seqs <- seqs[match(gt$tip.label, rownames(seqs)), ]
    slog("OK", file = logfile, megProj = x)
  }
  
  ## Write to database 
  ## -----------------
  slog("\nWriting alignment to database ... ", file = logfile, megProj = x)
  dbWriteMSA(x, dna = seqs, status = "aligned")
  dbDisconnect(conn)
  slog("OK", file = logfile, megProj = x)
  
  # Write files
  # -----------
  seqs <- dbReadMSA(x)
  slog("\nWriting alignment to PHYLIP file ... ", file = logfile, megProj = x)
  write.phy(seqs, paste0("msa/", gene, ".phy"))
  slog("OK", file = logfile, megProj = x)
  
  slog("\nWriting alignment to NEXUS file ... ", file = logfile, megProj = x)
  rownames(seqs) <- gsub("-", "_", rownames(seqs))
  write.nex(seqs, paste0("msa/", gene, ".nex"))
  slog("OK", file = logfile, megProj = x)

  ## calculate mean absolute deviation (MAD)
  ## see Smith, Beaulieau, Donoghue (2009)
  this.mad <- round(MAD(seqs), 5)
  slog("\nMean absolute deviation:",
       this.mad, "..", file = logfile, megProj = x)

  ## summary
  ## -------
  if (!quiet) {
    slog(
      paste("\n\n--- final alignment of", gene, "---"),
      paste("\nnumber of sequences     :", nrow(seqs)),
      paste("\nnumber of base pairs    :", ncol(seqs)),
      paste("\nmean absolute deviation :", this.mad),
      "\n\nSTEP G finished",
      file = logfile, megProj = x)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), "\n", file = logfile, megProj = x)
  }
  dbProgress(x, "step_g", "success")
}
