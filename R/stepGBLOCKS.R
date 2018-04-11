## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-01-29)

#' @title Detect Homology Uncertainty in Alingment
#' @description Use GBLOCKS to detect positions of uncertain homology and remove
#'   them from the alignment.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @return None, \code{stepGBLOCKS} is called for its side effects.
#' @seealso \code{\link{megapteraProj}}; \code{\link{stepMAFFT}} for the
#'   preceeding step and \code{\link{stepH}} for the subsequent step.
#' @export
#' @import DBI
#' @importFrom ips gblocks write.nex write.phy

stepGBLOCKS <- function(x){
  
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
    return()
  }
  ## Currently progress end with stepH (2018-02-26)
  # if (status$step_h == "success") {
  #   dbProgress(x, "step_i", "error")
  # }
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- x@taxon@tip.rank
  msa.tab <- paste(tip.rank, "sequence", sep = "_")
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
  logfile <- paste0("log/", gene, "-stepI.log")
  if (file.exists(logfile)) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP I: detecting nucleotide positions of uncertain homology\n",
       paste("\nLocus:", x@locus@sql),
       file = logfile)
  
  ## open database connection
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
  
  ## check if at least 3 species/genera are available
  ## -----------------------------------------
  n <- paste("SELECT count(taxon) FROM", msa.tab)
  n <- dbGetQuery(conn, n)$count
  if (n < 3){
    dbDisconnect(conn)
    slog("\nWARNING: only", n, "taxa available, no masking possible\n", 
         file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP I finished after", round(td, 2), attr(td, "units"), 
         "\n", file = logfile)
    return()
  }
  if (n < 100){ # 100 is arbitrary
    n <- paste("SELECT taxon FROM", msa.tab)
    n <- dbGetQuery(conn, n)
    n <- which(is.ingroup(x, n$taxon))
    if (length(n) < 3){
      dbDisconnect(conn)
      slog("\nWARNING:", length(n), "ingroup taxa available,", 
           "no masking possible\n", file = logfile)
      td <- Sys.time() - start
      slog("\nSTEP GG finished after", round(td, 2), attr(td, "units"), 
           "\n", file = logfile)
      return()
    }
  }
  
  ## read alignment
  ## --------------
  slog("\nReading alignment", file = logfile)
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
  a <- dbReadMSA(x)
  # }
  if (is.null(a)) {
    dbDisconnect(conn)
    slog("\nWARNING: no sequences conform to current parameter setting\n", file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP I finished after", round(td, 2), attr(td, "units"), "\n",
         file = logfile)
    return()
  }
  
  ## Masking of poorly aligned nucleotides
  ## -------------------------------------
  slog("\nMasking poorly aligned regions:", file = logfile)
  score <- gblocks(a, exec = gblocks.exe,
                   b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5, 
                   target = "score") # with least conservative default
  
  ## Write scores to database 
  ## ------------------------
  slog("\nWrite reliability scores to database ...", file = logfile)
  dbWriteMSA(x, dna = a, score = score, status = "aligned")
  slog("done", file = logfile)
  
  ## write files
  ## -----------
  # slog("\n.. write alignment to files ..", file = logfile)
  # if ( length(a) == 1 ) {
  #   write.phy(a[[1]], paste0("msa/", gene, "-masked.phy"))
  #   rownames(a[[1]]) <- gsub("-", "_", rownames(a[[1]]))
  #   write.nex(a[[1]], paste0("msa/", gene, "-masked.nex"))
  # } else {
  #   for ( i in 1:length(a) ){
  #     b <- paste("block", i, sep = "")
  #     write.phy(a[[i]], paste0("msa/", gene, "-", b, "-masked.phy"))
  #     rownames(a[[i]]) <- gsub("-", "_", rownames(a[[i]]))
  #     write.nex(a[[i]], paste0("msa/", gene, "-", b, "-masked.nex"))
  #   }
  #   a <- do.call(cbind.DNAbin, c(a, fill.with.gaps = TRUE))
  #   b <- paste(length(n), "blocks", sep = "")
  #   write.phy(a, paste0("msa/", gene, "-masked.phy"))
  #   rownames(a) <- gsub("-", "_", rownames(a))
  #   write.nex(a, paste0("msa/", gene, "-masked.nex"))
  # }
  
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
}
