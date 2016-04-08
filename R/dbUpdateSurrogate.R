## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-02-17)

dbUpdateSurrogate <- function(x, spec, surrogate){
  
  if ( !inherits(x, "megapteraProj") ) 
    stop("object 'x' is not of class 'megapteraProj'")
  if ( !is.character(spec) ) 
    stop("object 'spec' is of mode 'character'")
  if ( length(spec) == 2 ){
    surrogate <- spec[2]
    spec <- spec[1]
  }
  
  conn <- dbconnect(x)
  
  ## check if surrogate has been used already
  ## ----------------------------------------
  assigned <- paste("SELECT spec",
                    "FROM taxonomy",
                    "WHERE", wrapSQL(surrogate, term = "tag"))
  assigned <- dbGetQuery(conn, assigned)
  if ( nrow(assigned) > 0 ){
    if ( spec == assigned$spec ){
      dbDisconnect(conn)
      return()
    } else {
      stop(surrogate, " already assigned as surrogate for ", 
           assigned$spec)
    }
  }
  
  ## assign surrogate to spec
  ## ------------------------
  SQL <- paste("tag || ' surrogate:", surrogate, "'", sep = "")
  SQL <- paste("UPDATE taxonomy", 
               "SET tag =", SQL,
               "WHERE", wrapSQL(spec, operator = "="))
  dbSendQuery(conn, SQL)
  dbDisconnect(conn)
}