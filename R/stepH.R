## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2018-01-31)

#' @title Step H: Detect and Separate Unalignable Blocks
#' @description Dependent on the substitution rate of the genomic region and the
#'   taxonomic depth of the study group, not always all accession can be aligned
#'   into a single alignment. \code{step H} will detect and separate such
#'   unalignable alignment blocks.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param max.mad Numeric, giving the treshold value for the assessment of
#'   saturation: alignments with a median average distance (MAD) of
#'   \code{max.mad} or greater will be broken into blocks. The default value of
#'   0.01 has been estimated with simulation by Smith et al. (2009).
#' @return None, \code{stepH} is called for its side effects.
#' @references Smith, S. A., J. M. Beaulieu, and M. J. Donoghue. 2009.
#'   Mega-phylogeny approach for comparative biology: an alternative to
#'   supertree and supermatrix approaches. \emph{BMC Evolutionary Biology}
#'   \bold{9}:37.
#' @seealso \code{\link{megapteraProj}}; \code{\link{stepMAFFT}} for the preceeding
#'   step and \code{\link{supermatrix}} for the concatenation of loci.
#' @export
#' @import DBI snowfall

stepH <- function(x, max.mad){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("undefined locus not allowed")
  
  ## check if previous step has been run
  ## -----------------------------------
  status <- dbProgress(x)
  if (status$step_g == "pending") {
    stop("the previous step has not been called yet")
  }
  if (status$step_g == "error") {
    stop("the previous step has terminated with an error")
  }
  if (status$step_g == "failure") {
    slog("\nNo data from upstream available - quitting", file = "")
    dbProgress(x, "step_h", "failure")
    return()
  }
  if (status$step_g == "success") {
    dbProgress(x, "step_h", "error")
  }
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- x@taxon@tip.rank
  msa.tab <- paste(tip.rank, "sequence", sep = "_")
  block.max.dist <- x@params@block.max.dist
  min.n.seq <- x@params@min.n.seq
  if (missing(max.mad)) max.mad <- x@params@max.mad
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepH.log")
  if (file.exists(logfile)) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste0("\n", Sys.time()),
       "\nSTEP H: detection and separation of unsaturated blocks\n", 
       paste("\nLocus:", x@locus@sql), file = logfile)
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x)
  
  ## read alignment
  ## --------------
  slog("\nReading alignment with ", file = logfile)
  
  a <- dbReadMSA(x)
  
  ## check if stepG has been run properly
  ## should include checking for status = 'raw'
  if (!is.matrix(a)){
    stop("stepG has not been called yet")
  }
  ## how many species are in the table/alignment
  slog(nrow(a), "species", file = logfile)
  
  ## Set status column to 'aligned':
  ## If this is not done and the alignment is below
  ## MAD threshold, blocks from previous runs will
  ## not be eliminated.
  ## -----------------
  SQL <- paste("UPDATE", msa.tab,
               "SET status = 'aligned'",
               "WHERE status ~ 'block'",
               "AND", wrapSQL(gene, "locus", "="))
  dbSendQuery(conn, SQL)
  
  ## Checking saturation of alignment as measures
  ## by the median average difference (MAD) between
  ## corrected and JC69-corrected patristic distances
  ## ------------------------------------------------
  check.mad <- MAD(a)
  slog("\nAssessing saturation: MAD =", round(check.mad, 5),
       file = logfile)
  if (check.mad <= max.mad){
    SQL <- paste("UPDATE", msa.tab,
                 "SET status = '1 block'",
                 "WHERE status = 'aligned'",
                 "AND", wrapSQL(gene, "locus", "="))
    dbSendQuery(conn, SQL)
    dbDisconnect(conn)
    slog(paste("\nMAD below threshold of ", max.mad, ": ",
               "alignment not saturated\n", sep = ""), file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP H finished after", round(td, 2), attr(td, "units"), "\n",
         file = logfile)
    dbProgress(x, "step_h", "success")
    return()
  }
  
  ## Prepare guide tree 
  ## ------------------
  slog("\nPreparing guidetree ... ", file = logfile)
  gt <- comprehensiveGuidetree(x, tip.rank = tip.rank, subset = a)
  gt <- fixNodes(gt)
  a <- a[match(gt$tip.label, rownames(a)), ]
  a <- list(a)
  this.root <- Ntip(gt) + 1
  slog("OK", file = logfile)
  
  ## Split alignment into unsaturated blocks
  ## ---------------------------------------
  slog("\nSuccessively splitting up alignment:", file = logfile)
  repeat {
    aa <- splitAlignment(this.root, gt, a[[1]])
    a <- c(a[-1], aa)
    check.mad <- sapply(a, MAD)
    slog("\n- maximum MAD:", max(check.mad), file = logfile)
    if ( all(check.mad <= max.mad) ) break
    a <- a[order(check.mad, decreasing = TRUE)]
    this.root <- noi(gt, rownames(a[[1]]))
  }
  slog("\nAlignment was split into", length(a), 
       "unsaturated blocks", file = logfile)
  check.size <- sapply(a, nrow)
  a <- a[order(check.size, decreasing = TRUE)]
  
  ## Write results to database
  ## ----------------------------------
  tips <- lapply(a, function(z) gsub("_", " ", rownames(z)))
  mrca <- lapply(tips, taxdumpMRCA, x = x)
  names(tips) <- paste("block", 1:length(tips), mrca)
  tips <- lapply(tips, wrapSQL, term = "taxon", operator = "=", boolean = "OR")
  tips <- lapply(tips, function(z) paste0("(", z, ")"))
  SQL <- paste("UPDATE", msa.tab,
               "SET", wrapSQL(names(tips), "status", "=", NULL),
               "WHERE", tips,
               "AND", wrapSQL(gene, "locus", "="))
  lapply(SQL, dbSendQuery, conn = conn)
  
  ## summary of result
  ## -----------------
  blocks <- format(paste(sapply(a, nrow), "species"), justify = "right")
  blocks <- paste0("\n- ", format(names(tips)), blocks)
  slog(blocks, file = logfile)
  # slog(paste("\nMean absolute deviation :", this.mad),
  #      file = logfile)
  
  
  dbDisconnect(conn)
  slog("\n\nSTEP H finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), "\n", file = logfile)
  dbProgress(x, "step_h", "success")
}
