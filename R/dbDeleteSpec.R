## This code is part of the megaptera package
## © C. Heibl 2015 (last update: 2016-12-07)

#' @importFrom DBI dbDisconnect dbGetQuery dbSendQuery
#' @export

dbDeleteSpec <- function(x, spec){
  
  spec <- gsub(" ", "_", spec)
  
  ## taxonomy:
  conn <- dbconnect(x)
  SQL <- paste("DELETE FROM taxonomy",
                     "WHERE", wrapSQL(spec))
  dbSendQuery(conn, SQL)
  
  ## locus:
  SQL <- paste("DELETE FROM locus",
               "WHERE", wrapSQL(spec))
  dbSendQuery(conn, SQL)
  
  ## acc + spec:
  tabnames <- paste("SELECT table_name",
                    "FROM information_schema.tables",
                    "WHERE table_schema='public'",
                    "AND table_type='BASE TABLE'")
  tabnames <- dbGetQuery(conn, tabnames)$table_name
  acc <- tabnames[grep("^acc_", tabnames)]
  if ( length(acc) > 0 ){
    SQL <- paste("DELETE FROM", acc,
                 "WHERE", wrapSQL(spec, "taxon", operator = "="))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  spectab <- tabnames[grep("^spec", tabnames)]
  if ( length(spectab) > 0 ){
    SQL <- paste("DELETE FROM", spectab,
                 "WHERE", wrapSQL(spec, "spec", operator = "="))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  dbDisconnect(conn)
}