## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-04-11)

#' @export

dbUpdateReference_Clade <- function(megProj, tax){
  
  ## CHECKS
  ## ------
  if (!inherits(megProj, "megapteraProj"))
    stop("'megProj' is not of class 'megapteraProj'")
  if (megProj@locus@kind == "undefined") stop("undefined locus not allowed")
  if (!is.data.frame(tax)) stop("'tax' is not a data frame")
  names(tax) <- gsub("ref", "reference", names(tax))
  names(tax) <- gsub("tax$", "taxon", names(tax))
  if (!all(c("reference", "taxon") %in% names(tax)))
    stop("tax must contain columns 'reference' and 'taxon'")
  
  conn <- dbconnect(megProj)
  
  ## create table or delete entries for current locus
  ## ------------------------------------------------
  if (!dbExistsTable(conn, "reference_clade")) {
    SQL <- paste("CREATE TABLE reference_clade",
                 "(gene text NOT NULL,",
                 "taxon text NOT NULL,",
                 "reference text NOT NULL,",
                 "CONSTRAINT reference_clade_pk PRIMARY KEY (gene, taxon))")
    dbSendQuery(conn, SQL)
  } else {
    SQL <- paste("DELETE FROM reference_clade",
                 "WHERE", wrapSQL(megProj@locus@sql, "gene", "="))
    dbSendQuery(conn, SQL)
  }
  
  ## write data to table
  ## -------------------
  SQL <- paste(wrapSQL(megProj@locus@sql, NULL, "=", NULL),
               wrapSQL(tax$taxon, NULL, "=", NULL),
               wrapSQL(tax$reference, NULL, "=", NULL), sep = ", ")
  SQL <- paste0("(", SQL, ")")
  SQL <- paste("INSERT INTO reference_clade (gene, taxon, reference)",
               "VALUES", SQL)
  lapply(SQL, dbSendQuery, conn = conn)
  dbDisconnect(conn)
}