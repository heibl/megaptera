## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2019-02-03)

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
  write.fas(x@locus@reference, refn)
  cmd <- paste("/usr/local/ncbi/blast/bin/makeblastdb", 
               "-in", refn, 
               "-dbtype nucl",
               "-parse_seqids")
  system(cmd)
  
  ## Prepare Sequences
  ## ------------------
  seqs <- dbReadDNA(x, acc.tab)
  names(seqs) <- gsub("-", "__", names(seqs)) # hyphen illegal in BLAST!
  write.fas(seqs, dbn)
  
  ## Do the BLAST
  ## ------------
  cls <- c("qseqid", "sseqid", 
           "length", 
           "mismatch", 
           "qstart", "qend", "sstart", "send",
           # "gapopen", # Number of gap openings
           "qcovs", # Query Coverage Per Subject
           "pident", # Percentage of identical matches
           "evalue", 
           "bitscore", 
           "sstrand")
  outfmt <- paste0("-outfmt '", paste(c(6, cls), collapse = " "), "'")
  cmd <- paste("/usr/local/ncbi/blast/bin/blastn",
               "-db", refn,
               "-query", dbn,
               "-task blastn", ## traditional BLASTN requiring an exact match of 11
               "-evalue 1000",
               "-max_hsps 1", ## maximum number of HSPs per subject sequence to save for each query
               "-max_target_seqs 10000",
               outfmt,
               "-out", outn)
  system(cmd)
  
  ## Parse output 
  ## ------------
  res <- fread(outn)
  names(res) <- cls
  ## For every species select hit with lowest E-value
  res <- by(res, res$qseqid, function(z) z[which.min(z$evalue),])
  res <- do.call(rbind, res)
  res$qseqid <- gsub("__", "-", res$qseqid) ## bring hyphen back!
  res <- cbind(splitGiTaxon(res$qseqid), res)
  
  # test <- res[grep("Miniopterus_pallidus", res$qseqid), ]
  # test <- res[grep("Balae", res$qseqid), ]
  
  ## Update database
  ## ---------------
  SQL <- paste("UPDATE", acc.tab,
               "SET", wrapSQL(round(res$evalue, 10), "e_value", "=", NULL),
                ",", wrapSQL(res$pident, "identity", "=", NULL),
                ",", wrapSQL(res$qcovs, "coverage", "=", NULL),
               "WHERE", wrapSQL(res$gi, "gi", "=", NULL))
  conn <- dbconnect(x)
  lapply(SQL, dbSendQuery, conn = conn)
  dbDisconnect(conn)
  
}