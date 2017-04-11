## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-10-28)

#' @export

dbconnect <- function(dbPars){
  
  if ( !inherits(dbPars, c("dbPars", "megapteraProj")) ){
    stop("dbPars is not of classes 'dbPars' or 'megapteraProj'")
  }
  if ( inherits(dbPars, "megapteraProj") ){
    dbPars <- dbPars@db
  }
  DBI::dbConnect(RPostgreSQL::PostgreSQL(), 
                 host = dbPars@host,
                 port = dbPars@port,
                 user = dbPars@user, 
                 password = dbPars@password, 
                 dbname = dbPars@dbname)
}