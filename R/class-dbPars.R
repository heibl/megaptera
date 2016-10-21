## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-10-21)

setClass("dbPars", 
         representation = list(
           host = "character", 
           port = "numeric", 
           dbname = "character", 
           user = "character", 
           password = "character")
)

## USER LEVEL CONSTRUCTOR
## ----------------------
"dbPars" <- function(host = "localhost", 
                     port = 5432, 
                     dbname, 
                     user = "postgres", 
                     password){
  
  dbname <- tolower(dbname)
  
  ## check if database exists ...
  conn <- DBI::dbConnect(RPostgreSQL::PostgreSQL(), 
                    dbname = "postgres",
                    host = host,
                    user = user, 
                    port = port, 
                    password = password)
  sql <- paste("SELECT 1 FROM pg_database WHERE",
               sql.wrap(dbname, term = "datname"))
  if ( nrow(dbGetQuery(conn, sql)) == 1 ){
    cat("\ndatabase '", dbname, "' exists", sep = "")  
  } else {
    ## .. and create if it does not exist
    cat("\ndatabase '", dbname, "' created", sep = "") 
    sql <- paste("CREATE DATABASE", dbname,
                 "WITH ENCODING='UTF8'",
                 "CONNECTION LIMIT=-1;")
    dbSendQuery(conn, sql)
  }
  dbDisconnect(conn)
  
  new("dbPars", 
      host = host, port = port, 
      dbname = dbname, 
      user = user, password = password
  )
}

## SET METHOD: SHOW
## ----------------
setMethod("show",
          signature(object = "dbPars"),
          function (object) 
          {
            cat("PostgreSQL connection parameters:",
                "\n     host =", object@host,
                "\n     port =", object@port,
                "\n   dbname =", object@dbname,
                "\n     user =", object@user,
                "\n password =", object@password
            )
          }
)

## FUNCTION: DBCONNECT
## -------------------
dbconnect <- function(dbPars){
  if ( !inherits(dbPars, c("dbPars", "megapteraProj")) ){
    stop("dbPars is not of classes 'dbPars' or 'megapteraProj'")
  }
  if ( inherits(dbPars, "megapteraProj") ){
    dbPars <- dbPars@db
  }
  dbConnect(PostgreSQL(), 
            host = dbPars@host,
            port = dbPars@port,
            user = dbPars@user, 
            password = dbPars@password, 
            dbname = dbPars@dbname)
}