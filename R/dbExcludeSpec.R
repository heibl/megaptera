## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2016-12-07)

#' @importFrom DBI dbDisconnect dbSendQuery
#' @export

dbExcludeSpec <- function(x, spec){
  
  gene <- x@locus@sql
  tip.rank <- x@taxon@tip.rank
  msa.tab <- paste(tip.rank, gsub("^_", "", gene), sep = "_")

  ## update spec_* table
  ## -------------------
  conn <- dbconnect(x)
  SQL <- paste("UPDATE", msa.tab,
               "SET status='excluded (by user)'",
               "WHERE", wrapSQL(spec, 
                                term = tip.rank, 
                                operator = "=", 
                                boolean = NULL))
  lapply(SQL, dbSendQuery, conn = conn)
  SQL <- paste("UPDATE", msa.tab,
               "SET status='raw'",
               "WHERE status!='excluded (by user)'")
  dbSendQuery(conn, SQL)
  dbDisconnect(conn)
}