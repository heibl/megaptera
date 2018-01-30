## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2018-01-29)

#' @export
#' @import DBI

dbWriteMSA <- function(megProj, dna, 
                       status = "raw",
                       n, md5, score = FALSE){
  
  gene <- megProj@locus@sql
  tip.rank <-megProj@taxon@tip.rank
  msa.tab <- paste(tip.rank, "sequence", sep = "_")

  conn <- dbconnect(megProj)
  
  ## Prepare data
  ## ------------
  dna <- DNAbin2pg(dna, score)
  dna <- cbind(locus = gene, dna, status = status)
  if (!missing(n)) dna <- cbind(dna, n)
  if (!missing(md5)) dna <- cbind(dna, md5)
  dna$taxon <- gsub("_", " ", dna$taxon)
  
  ## Devide into present and new taxa
  ## ---------------------------------
  taxa_present <- paste("SELECT taxon FROM", msa.tab,
                        "WHERE", wrapSQL(gene, "locus", "="))
  taxa_present <- dbGetQuery(conn, taxa_present)$taxon
  taxa_present <- intersect(dna$taxon, taxa_present)
  taxa_new <- setdiff(dna$taxon, taxa_present)
  
  ## STORE DATA:
  
  ## 1. Add meta data for new taxa
  ## -----------------------------
  # cat("\nAdding", length(taxa_new), "new taxa")
  if (length(taxa_new)){
    dna_new <- dna[dna$taxon %in% taxa_new, ]
    SQL <- apply(dna_new, 1, function(z) paste(paste0("'", z, "'"), collapse = ", "))
    SQL <- paste(paste0("(", SQL, ")"), collapse = ", ")
    SQL <- paste("INSERT INTO", msa.tab, 
                 paste0("(", paste(names(dna_new), collapse = ", "), ")"),  
                 "VALUES", SQL)
    dbSendQuery(conn, SQL)
  }
  
  ## 2. Update meta data and delete sequence 
  ##    data for present taxa
  ## ---------------------------------------
  # cat("\nUpdating", length(taxa_present), "present taxa")
  if (length(taxa_present)){
    tt <- dna[dna$taxon %in% taxa_present, ]
    SQL <- tt[, !names(tt) %in% c("locus", "taxon")]
    for (i in 1:ncol(SQL)){
      SQL[, i] <- paste0(names(SQL)[i], "='", SQL[, i], "'")
    }
    SQL <- apply(SQL, 1, paste, collapse = ", ")
    
    SQL <- paste("UPDATE", msa.tab, 
                 "SET", SQL, 
                 "WHERE", wrapSQL(gene, "locus", "="),
                 "AND", wrapSQL(tt$taxon, "taxon", "=", NULL))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  
  ## Close database connection
  ## -------------------------
  invisible(dbDisconnect(conn))
}
