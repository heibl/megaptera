## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2016-12-07)

#' @importFrom DBI dbDisconnect dbSendQuery
#' @export

dbIncludeSpec <- function(x, spec){
  
  gene <- x@locus@sql
  tip.rank <- x@taxon@tip.rank
  msa.tab <- paste(tip.rank, gsub("^_", "", gene), sep = "_")

  ## update spec_* table
  ## -------------------
  conn <- dbconnect(x)
  SQL <- paste("UPDATE", msa.tab,
               "SET status='raw'")
  dbSendQuery(conn, SQL)
  dbDisconnect(conn)
}