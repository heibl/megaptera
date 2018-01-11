## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-02-17)

#' @export
#' @import DBI

stepPASTAin <- function(x){	
  
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
  
  ## open database connection
  conn <- dbconnect(x)
  
  ## check if msa table exists
  ## -------------------------
  if ( !dbExistsTable(conn, msa.tab) ){
    dbDisconnect(conn)
    slog("\nWARNING: table", msa.tab, "does not exist!\n")
    td <- Sys.time() - start
    slog("\nSTEP G finished after", round(td, 2), attr(td, "units"), 
         "\n")
    return()
  }
  dbDisconnect(conn)
  
  ## sequences
  ## ---------
  seqs <- dbReadDNA(x, msa.tab)
  seqs <- del.gaps(seqs)
  
  ## prepare guide tree for alignment
  ## --------------------------------
  gt <- comprehensiveGuidetree(x, tip.rank = "spec", subset = seqs)
  
  ## name checking
  ## -------------
  if (!all(gt$tip.label %in% names(seqs)) | !all(names(seqs) %in% gt$tip.label))
 
  ## write files
  ## -----------
  write.fas(seqs, paste(gene, "fas", sep = "."))
  write.tree(gt, paste(gene, "tre", sep = "."))
  
  cat("cd", getwd())
  cat("\npython", 
      "/Applications/pasta-code/pasta/run_pasta.py",
      "-i", paste(gene, "fas", sep = "."),
      "-t", paste(gene, "tre", sep = "."),
      "-j", gene)
}
