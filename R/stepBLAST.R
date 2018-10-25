## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2018-04-12)

#' @title stepBLASTN: Filter Homologous Sequences
#' @description Filters homologous sequences using the BLASTB tool
#' @param x An object of class \code{\link{megapteraProj}}.
#' @references NCBI BLAST: \url{https://blast.ncbi.nlm.nih.gov/Blast.cgi}
#' @importFrom data.table fread
#' @export

stepBLAST <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("undefined locus not allowed")
  
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  dbn <- paste0("data/acc_", gsub("^_", "", gene), ".fas")
  refn <- paste0("data/ref_", gsub("^_", "", gene), ".fas")
  outn <- paste0("data/out_", gsub("^_", "", gene), ".txt")
  # dbProgress(x, "step_b", "error")
  
  ## Create BLAST database
  ## ---------------------
  seqs <- dbReadDNA(x, acc.tab)
  write.fas(seqs, dbn)
  cmd <- paste("/usr/local/ncbi/blast/bin/makeblastdb", 
               "-in", dbn, 
               "-dbtype nucl",
               "-parse_seqids")
  system(cmd)
  
  ## Prepare References
  ## ------------------
  write.fas(x@locus@reference, refn)
  
  ## Do the BLAST
  ## ------------
  cls <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
           "gapopen", "qstart", "qend", "sstart", "send", 
           "evalue", "bitscore")
  outfmt <- paste0("-outfmt '", paste(c(6, cls), collapse = " "), "'")
  cmd <- paste("/usr/local/ncbi/blast/bin/blastn",
               "-db", dbn,
               "-query", refn,
               "-task blastn-short",
               "-evalue 1000",
               "-max_target_seqs 10000",
               outfmt,
               "-out", outn)
  system(cmd)
  
  ## Parse output and update database
  ## --------------------------------
  res <- fread(outn)
  names(res) <- cls
  res <- cbind(splitGiTaxon(res$sseqid), res)
  
  SQL <- paste("UPDATE", acc.tab,
               "SET", wrapSQL(round(res$evalue, 10), "e_value", "=", NULL),
                ",", wrapSQL(res$pident, "identity", "=", NULL),
                ",", wrapSQL(res$length, "coverage", "=", NULL),
               "WHERE", wrapSQL(res$gi, "gi", "=", NULL))
  conn <- dbconnect(x)
  lapply(SQL, dbSendQuery, conn = conn)
  dbDisconnect(conn)
  
}