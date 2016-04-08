## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2014-12-19)

icContingency <- function(x, step = .05){
  
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  min.identity <- x@locus@min.identity
  min.coverage <- x@locus@min.coverage
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
  
  SQL <- paste("SELECT taxon,",
               "max(identity) AS maxi, max(coverage) AS maxc",
               "FROM", acc.tab, 
               "WHERE identity IS NOT NULL",
               "AND coverage IS NOT NULL",
               "AND status !~ 'excluded|too'",
               "GROUP BY taxon")
  z <- dbGetQuery(conn, SQL)
  dbDisconnect(conn)
  
  mi <- seq(0, 1, step)
  mc <- seq(0, 1, step)
  obj <- matrix(nrow = length(mi), ncol = length(mc), 
                dimnames = list(mi, mc))
  
  for ( i in seq_along(mi)){
    for ( j in seq_along(mc)){
      obj[i, j] <- nrow(z[z$maxi >= mi[i] & z$maxc >= mc[j], ])
    }
  }
  obj
}