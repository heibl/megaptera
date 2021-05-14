 ## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2021-03-19)

#' @export
#' @import DBI

dbWriteMSA <- function(megProj, dna, 
                       status = "raw",
                       n, md5, score = FALSE){
  
  locus <- megProj@locus@sql
  if (megProj@locus@kind == "undefined") {
    stop("locus undefined; use setLocus() to define a locus")
    }
  tip.rank <- megProj@taxon@tip.rank
  msa.tab <- "sequence_selected"

  conn <- dbconnect(megProj)
  
  ## Prepare data
  ## ------------
  dna <- DNAbin2pg(dna, score)
  dna <- cbind(locus = locus, dna, status = status)
  if (!missing(n)) dna <- cbind(dna, n)
  if (!missing(md5)) dna <- cbind(dna, md5)
  label <- splitGiTaxon(dna$taxon, white.space = " ")
  dna$taxon <- label$taxon
  if (!is.null(label$gi)) dna$acc <- label$gi
  
  ## Divide into present and new taxa
  ## ---------------------------------
  taxa_present <- paste("SELECT taxon FROM sequence_selected", ## ||' '||locus 
                        "WHERE", wrapSQL(locus, "locus", "="))
  taxa_present <- dbGetQuery(conn, taxa_present)[, 1]
  if (any(duplicated(taxa_present))) stop("debug me! [dbWriteMSA l. 35]")
  dna$new <- !dna$taxon %in% taxa_present
  
  ## STORE DATA:
  
  ## 1. Add meta data for new taxa
  ## -----------------------------
  # cat("\nAdding", length(taxa_new), "new taxa")
  if (any(dna$new)){
    dna_new <- dna[dna$new, names(dna) != "new"]
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
  if (!all(dna$new)){
    tt <- dna[!dna$new, names(dna) != "new"]
    SQL <- tt[, !names(tt) %in% c("locus", "taxon", "acc", "new")]
    for (i in 1:ncol(SQL)){
      SQL[, i] <- paste0(names(SQL)[i], "='", SQL[, i], "'")
    }
    SQL <- apply(SQL, 1, paste, collapse = ", ")
    
    SQL <- paste("UPDATE", msa.tab, 
                 "SET", SQL, 
                 "WHERE", wrapSQL(locus, "locus", "="),
                 "AND", wrapSQL(tt$taxon, "taxon", "=", NULL))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  
  ## Close database connection
  ## -------------------------
  invisible(dbDisconnect(conn))
}
