## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2020-02-20)

#' @title stepBLASTN: Filter Homologous Sequences
#' @description Filters homologous sequences using the BLASTB tool
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param word.size Integer, ...
#' @param gapopen Integer, ...
#' @param gapextend Integer, ...
#' @param penalty Integer, ...
#' @param reward Integer, ...
#' @details 
#' \tabular{ll}{
#' \code{sseqid} \tab ID of query sequence that gave the lowest E-value \cr
#' \code{length} \tab Alignment length \cr
#' \code{mismatch} \tab Alignment length \cr
#' \code{qstart} \tab Start of alignment in query \cr
#' \code{qend} \tab End of alignment in query \cr
#' \code{sstart} \tab Start of alignment in subject \cr
#' \code{send} \tab End of alignment in subject \cr
#' \code{qcovs} \tab Percentage of query coverage per subject\cr
#' \code{pident} \tab Percentage of identical matches \cr
#' \code{evalue} \tab Expect value \cr
#' \code{bitscore} \tab Bit score \cr
#' \code{sstrand} \tab Percentage of identical matches \cr
#' }
#' @references NCBI BLAST: \url{https://blast.ncbi.nlm.nih.gov/Blast.cgi}
#' @importFrom data.table fread
#' @export

stepBLAST <- function(x, word.size = 11, gapopen = 5, gapextend = 2,
                      penalty = -3, reward = 1){
  
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
  conn <- dbconnect(x)
  seqs <- dbGetQuery(conn, "SELECT acc, sequence FROM sequence")
  acc_table <- data.frame(raw = seqs$acc,
                          canonical = gsub(" |-", "_", seqs$acc),
                          stringsAsFactors = FALSE) ## space and hyphen not allowed!
  seqs$sequence <- gsub("-", "", seqs$sequence) ## delete gap characters
  seqs <- paste(paste0(">", acc_table$canonical), seqs$sequence, sep = "\n")
  cat(seqs, file = dbn, sep = "\n")
  remove(seqs)
  
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
               "-task blastn",  ## traditional BLASTN requiring an exact match of 11
               "-evalue 10",
               "-word_size", word.size,   ## default: 11
               "-gapopen", gapopen,       ## default: 5
               "-gapextend", gapextend,   ## default: 2
               "-penalty", penalty,       ## default: -3
               "-reward", reward,         ## default: 1
               "-max_hsps 1",   ## maximum number of HSPs per subject sequence to save for each query
               "-max_target_seqs 10000",
               outfmt,
               "-out", outn)
  system(cmd)
  
  ## Parse output 
  ## ------------
  res <- fread(outn)
  names(res) <- cls
  ## For every accession select hit with lowest E-value
  res <- by(res, res$qseqid, function(z) z[which.min(z$evalue),])
  res <- do.call(rbind, res)
  
  ## Revert to raw accession numbers
  res$qseqid <- acc_table$raw[match(res$qseqid, acc_table$canonical)]
  
  message("blastn returned ", nrow(res), " hits (E-value <=", max(res$evalue), ")")
  
  
  ## We only keep hits which improve previous E-values
  ## -------------------------------------------------
  previous <- "SELECT acc, qcovs, evalue FROM sequence WHERE locus IS NOT NULL"
  previous <- dbGetQuery(conn, previous)
  previous <- previous[previous$acc %in% res$qseqid, ]
  
  res <- data.frame(res, previous[match(res$qseqid, previous$acc), ])
  res$evalue.1[is.na(res$evalue.1)] <- 1000 ## these are the accs not yet BLASTED
  res <- res[res$evalue < res$evalue.1, -(14:16)]
  message(nrow(res), " of these hits are new or represent lower E-values")
  
  ## Update database
  ## ---------------
  if (nrow(res)){
    SQL <- paste("UPDATE sequence",
                 "SET", wrapSQL(gene, "locus", "=", NULL),
                 ",", wrapSQL(res$sseqid, "sseqid", "=", NULL),
                 ",", wrapSQL(res$length, "length", "=", NULL),
                 ",", wrapSQL(res$mismatch, "mismatch", "=", NULL),
                 ",", wrapSQL(res$qstart, "qstart", "=", NULL),
                 ",", wrapSQL(res$qend, "qend", "=", NULL),
                 ",", wrapSQL(res$sstart, "send", "=", NULL),
                 ",", wrapSQL(res$qcovs, "qcovs", "=", NULL), ## coverage
                 ",", wrapSQL(res$pident, "pident", "=", NULL), ## identity
                 ",", wrapSQL(round(res$evalue, 10), "evalue", "=", NULL),
                 ",", wrapSQL(res$bitscore, "bitscore", "=", NULL),
                 ",", wrapSQL(res$sstrand, "sstrand", "=", NULL),
                 "WHERE", wrapSQL(res$qseqid, "acc", "=", NULL))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  dbDisconnect(conn)
}
