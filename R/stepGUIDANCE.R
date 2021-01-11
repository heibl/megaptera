## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-11-07)

#' @title Detect Homology Uncertainty in Alingment
#' @description Use GUIDANCE, GUIDANCE2 or HoT to calculate column-wise
#'   reliability scores for alignments.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param bootstrap An integer giving the number of alternative MSAs to be computed.
#' @return None, \code{stepGUIDANCE} is called for its side effects.
#' @seealso \code{\link{megapteraProj}}; \code{\link{stepMAFFT}} for the
#'   preceeding step and \code{\link{stepH}} for the subsequent step.
#' @export
#' @import DBI
#' @importFrom ips gblocks write.nex write.phy
#' @importFrom rGUIDANCE guidance scores

stepGUIDANCE <- function(x, bootstrap = 100){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("undefined locus not allowed")
  
  ## check if previous step has been run
  ## -----------------------------------
  # status <- dbProgress(x)
  # if (status$step_g == "pending") {
  #   stop("the previous step has not been called yet")
  # }
  # if (status$step_g == "error") {
  #   stop("the previous step has terminated with an error")
  # }
  # if (status$step_g == "failure") {
  #   slog("\nNo data from upstream available - quitting", file = "")
  #   return()
  # }
  ## Currently progress end with stepH (2018-02-26)
  # if (status$step_h == "success") {
  #   dbProgress(x, "step_i", "error")
  # }
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- x@taxon@tip.rank
  msa.tab <- "sequence_selected"
  block.max.dist <- x@params@block.max.dist
  max.bp <- x@params@max.bp * 1.5
  min.n.seq <- x@params@min.n.seq
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepI.log")
  if (file.exists(logfile)) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP GUIDANCE: assesssing reliability of nucleotide positions of uncertain homology\n",
       paste("\nLocus:", x@locus@sql),
       file = logfile)
  
  ## Open database connection
  ## ------------------------
  conn <- dbconnect(x)
  
  ## check if msa table exists
  ## -------------------------
  if (!dbExistsTable(conn, msa.tab)){
    dbDisconnect(conn)
    slog("\nWARNING: table", msa.tab, "does not exist!\n", file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP I finished after", round(td, 2), attr(td, "units"), "\n",
         file = logfile)
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
  
  ## Read alignment
  ## --------------
  slog("\nReading alignment", file = logfile)
  a <- dbReadMSA(x, confid.scores = "ignore")
  if (is.null(a)) {
    dbDisconnect(conn)
    slog("\nWARNING: no sequences conform to current parameter setting\n", file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP I finished after", round(td, 2), attr(td, "units"), "\n",
         file = logfile)
    return()
  }
  
  ## Calculating GUIDANCE residue scores
  ## Note: Column scores that are calculated from residue
  ## scores are mathematical identical to GUIDANCE CS
  ## ------------------------------------------------
  slog("\nCalculating column reliability score:\n", file = logfile)
  s <- guidance(a, msa.exec = x@align.exe, ncore = x@params@cpus, 
                bootstrap = bootstrap)
  a <- s@msa
  s <- scores(s, score = "residue", na.rm = FALSE)$residue
  slog("done", file = logfile)
  if (!identical(dim(a), dim(s))) stop("MSA and scores do not match; this is likely a bug")
  slog("\nMSA with", nrow(a), "rows and", ncol(a), "columns",  file = logfile)
  id <- is.na(s)
  s[id] <- 0 ## could also be the column average
  usn <- table(id)["TRUE"]
  slog("\n", usn, " unscored nucleotides (", round(usn / prod(dim(id)) * 100, 2),
       "%) are set to GUIDANCE residue score = 0", sep = "", file = logfile)
  slog("\nRange of scores:", paste(range(s), collapse = " - "), file = logfile)
  
  ## Write scores to database 
  ## ------------------------
  slog("\nWrite reliability scores to database ...", file = logfile)
  dbWriteMSA(x, dna = a, score = s, status = "aligned + reliability")
  slog("done", file = logfile)
  
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
  
  
  dbDisconnect(conn)
  slog("\n\nSTEP I finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
}
