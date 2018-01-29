## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2018-01-29)

#' @export
#' @import DBI

dbWriteMSA <- function(megProj, dna, 
                       status = "raw",
                       n, md5, score = FALSE){
  
  gene <- megProj@locus@sql
  tip.rank <- megProj@taxon@tip.rank
  msa.tab <- "species_meta"
  msa.seq.tab <- "species_sequence"
  nmd5 <- ifelse(missing(md5), "",
                 paste(wrapSQL(n, "n", "="), ",", 
                       wrapSQL(md5, "md5", "="), ","))
  
  conn <- dbconnect(megProj)
  
  ## Prepare data
  ## ------------
  dna <- DNAbin2pg(dna)
  taxa <- unique(dna$taxon)
  
  ## Devide into present and new taxa
  ## ---------------------------------
  taxa_present <- dbGetQuery(conn, paste("SELECT taxon FROM", msa.tab))$taxon
  taxa_new <- setdiff(taxa, taxa_present)
  
  ## STORE DATA:
  
  ## 1. Add meta data for new taxa
  ## -----------------------------
  if (length(taxa_new)){
    SQL <- paste("INSERT INTO", msa.tab, 
                 "(locus, taxon, n, md5, status)",  
                 "VALUES (", 
                 wrapSQL(gene, NULL), ",",
                 wrapSQL(taxa_new, NULL), ",",
                 wrapSQL(n, NULL), ",",
                 wrapSQL(md5, NULL), ",",
                 wrapSQL(status, NULL), ")")
    dbSendQuery(conn, SQL)
  }
  
  ## 2. Update meta data and delete sequence 
  ##    data for present taxa
  ## ---------------------------------------
  if (length(taxa_present)){
    SQL <- paste("UPDATE", msa.tab, 
                 "SET", 
                 nmd5, 
                 wrapSQL(status, "status", "="), ",", 
                 "WHERE", wrapSQL(gene, "locus"),
                 "AND", wrapSQL(taxa_present, "taxon"))
    dbSendQuery(conn, SQL)
    
    SQL <- paste("DELETE FROM", msa.seq.tab, 
                 "WHERE", wrapSQL(gene, "locus"),
                 "AND", wrapSQL(taxa_present, "taxon"))
    dbSendQuery(conn, SQL)
  }
  
  ## 3. Add sequence data
  ## --------------------
  SQL <- paste("INSERT INTO", msa.seq.tab, 
        "(locus, taxon, nuc, pos, reliability)",  
        "VALUES (", 
        wrapSQL(gene, NULL), ",",
        wrapSQL(dna$taxon, NULL, "=", NULL), ",",
        wrapSQL(dna$nuc, NULL, "=", NULL), ",",
        wrapSQL(dna$pos, NULL, "=", NULL), ",",
        wrapSQL(dna$reliability, NULL, "=", NULL), ")")
  lapply(SQL, dbSendQuery, conn = conn)
  
  ## Close database connection
  ## -------------------------
  dbDisconnect(conn)
}
