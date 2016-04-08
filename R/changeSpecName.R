## This code is part of the megaptera package
## © C. Heibl 2015 (last update: 2015-12-23)

dbChangeSpecName <- function(x, from, to){
 
  from <- gsub(" ", "_", from)
  to <- gsub(" ", "_", to)
  
  conn <- dbconnect(x)
  to.exists <- paste("SELECT synonym",
                     "FROM taxonomy",
                     "WHERE", wrapSQL(to))
  to.exists <- dbGetQuery(conn, to.exists)
  
  from.exists <- paste("SELECT synonym",
               "FROM taxonomy",
               "WHERE", wrapSQL(from))
  from.exists <- dbGetQuery(conn, from.exists)
  

  
  if ( nrow(to.exists) == 0 ){
    if ( from.exists$synonym == "-" ) synonym <- from
    else stop("synonym not empty - implement me!")
    tax <- paste("UPDATE taxonomy",
                 "SET", wrapSQL(strip.spec(to), "gen", operator = "="),
                 ",", wrapSQL(to, "spec", operator = "="),
                 ",", wrapSQL(synonym, "synonym", operator = "="),
                 "WHERE", wrapSQL(from, "spec", operator = "="))
    dbSendQuery(conn, tax)
  }
  
  ## acc + spec
  tabnames <- paste("SELECT table_name",
                    "FROM information_schema.tables",
                    "WHERE table_schema='public'",
                    "AND table_type='BASE TABLE'")
  tabnames <- dbGetQuery(conn, tabnames)$table_name
  acc <- tabnames[grep("^acc_", tabnames)]
  if ( length(acc) > 0 ){
    SQL <- paste("UPDATE", acc, 
                 "SET", wrapSQL(to, "taxon", operator = "="),
                 "WHERE", wrapSQL(from, "taxon", operator = "="))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  spec <- tabnames[grep("^spec", tabnames)]
  if ( length(spec) > 0 ){
    SQL <- paste("UPDATE", spec, 
                 "SET", wrapSQL(to, "spec", operator = "="),
                 "WHERE", wrapSQL(from, "spec", operator = "="))
    lapply(SQL, dbSendQuery, conn = conn)
  }
}