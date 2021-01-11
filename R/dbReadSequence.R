
#' @export
dbReadSequence <- function(x, locus = x@locus){
  
  conn <- dbconnect(x)
  seqs <- paste("SELECT taxon, acc, sequence", 
                "FROM sequence",
                "WHERE", wrapSQL(locus@sql, "locus", "="),
                "ORDER BY taxon")
  seqs <- dbGetQuery(conn, seqs)
  dbDisconnect(conn)
  seq_names <- gsub(" ", "_", paste(seqs$taxon, seqs$acc))
  seqs <- strsplit(seqs$sequence, "")
  names(seqs) <- seq_names
  as.DNAbin(seqs)
}

#' @export
dbReadSequenceSelected <- function(x, locus = x@locus@sql, status = "reliability", threshold = 0.25){
  
  ## Query database
  ## --------------
  conn <- dbconnect(x)
  seqs <- paste("SELECT taxon, sequence, reliability", 
                "FROM sequence_selected",
                "WHERE", wrapSQL(locus, "locus", "="),
                "AND", wrapSQL("reliability", "status", "~"),
                "ORDER BY taxon")
  seqs <- dbGetQuery(conn, seqs)
  dbDisconnect(conn)
  
  ## Reorganize data
  ## ---------------
  seq_names <- gsub(" ", "_", paste(seqs$taxon))
  reliability_scores <- strsplit(seqs$reliability, " ")
  reliability_scores <- lapply(reliability_scores, as.numeric)
  names(reliability_scores) <- seq_names
  seqs <- strsplit(seqs$sequence, "")
  names(seqs) <- seq_names
  
  ## Set zero weigth characters to 'N'
  ## ---------------------------------
  for (i in seq_along(reliability_scores)){
    zero <- reliability_scores[[i]] < threshold
    seqs[[i]][zero] <- "n"
  }
  
  ## Convert to DNAbin
  ## -----------------
  seqs <- as.DNAbin(seqs)
  if (nrow(unique(sapply(seqs, nchar))) == 1){
    seqs <- as.matrix(seqs)
  }
  # attr(seqs, "reliability_scores") <- reliability_scores
  seqs
}