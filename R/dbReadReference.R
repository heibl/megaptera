## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-11-03)

#' @export
#' @import DBI

dbReadReference <- function(x, locus){
  
  conn <- dbconnect(x)
  if (missing(locus)) locus <- x@locus@sql
  
  if (!dbExistsTable(conn, "reference")){
    dbDisconnect(conn)
    return(FALSE)
  }
  sql <- paste("SELECT taxon, reference", 
               "FROM reference",
               "WHERE", wrapSQL(locus, "gene", "="))
  sql <- dbGetQuery(conn, sql)
  dbDisconnect(conn)
  if ( nrow(sql) > 0 ){
    b <- as.DNAbin(strsplit(sql$reference, ""))
    names(b) <- sql$taxon
    if ( length(unique(sapply(b, length))) == 1 )
      b <- as.matrix(b)
  } else {
    b <- FALSE
  }
  b
}