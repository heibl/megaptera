## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-01-11)

#' @title Set Database Parameters
#' @description Sets the connection parameters for a PostgreSQL database and, as
#'   a side effect, creates the database, if it does not yet exist.
#' @param host A vector of mode \code{"character"}, giving the name of the host,
#'   most often this will be \code{"localhost"}.
#' @param port Numeric, giving the port number, most often \code{5432}.
#' @param dbname A vector of mode \code{"character"}, giving the name of the
#'   database.
#' @param user A vector of mode \code{"character"}, giving the name of the user.
#' @param password A vector of mode \code{"character"}, giving the password.
#' @details \bold{megaptera} stores all data internally in a PostgreSQL
#'   database. Therefore, \code{dbPars} represents the first step in setting up a megaptera
#'   project pipeline. See \code{\link[DBI]{dbConnect}} and
#'   \code{\link{PostgreSQL}} for further details about the connection
#'   parameters and procedure.
#' @references See the documentation at the PostgreSQL Web site
#'   \url{http://www.postgresql.org} for details.
#' @return An object of class \code{\linkS4class{dbPars}}.
#' @seealso \code{\linkS4class{dbPars}} for the class' description;
#'   \code{\link{taxon}}, \code{\link{locus}}, and \code{\link{megapteraPars}}
#'   for defining of taxa, loci, and the pipeline's parameters, respectively;
#'   and \code{\link{megapteraProj}} for the bundling of input data.
#' @include dbPars-class.R
#' @importFrom methods new
#' @export

## USER LEVEL CONSTRUCTOR
## ----------------------
"dbPars" <- function(host = "localhost", 
                     port = 5432, 
                     dbname, 
                     user = "postgres", 
                     password){
  
  dbname <- tolower(dbname)
  
  ## check if database exists ...
  ## ----------------------------
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
    ## ----------------------------------
    cat("\ndatabase '", dbname, "' created", sep = "") 
    sql <- paste("CREATE DATABASE", dbname,
                 "WITH ENCODING='UTF8'",
                 "CONNECTION LIMIT=-1;")
    dbSendQuery(conn, sql)
  }
  dbDisconnect(conn)
  
  ## create relation 'progress' if necessary
  ## ---------------------------------------
  conn <- DBI::dbConnect(RPostgreSQL::PostgreSQL(), 
                         dbname = dbname,
                         host = host,
                         user = user, 
                         port = port, 
                         password = password)
  if (!dbExistsTable(conn, "progress")){
    SQL <- paste("CREATE TABLE progress",
                 "(",
                 "locus text NOT NULL,",
                 "step_b text,",
                 "step_c text,",
                 "step_d text,",
                 "step_e text,",
                 "step_f text,",
                 "step_g text,",
                 "step_h text,",
                 "step_i text,",
                 "CONSTRAINT progress_pk PRIMARY KEY (locus)",
                 ")")
    dbSendQuery(conn, SQL)
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
