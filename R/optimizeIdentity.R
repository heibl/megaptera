## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2016-12-06)

#' @importFrom graphics hist
#' @export

optimizeIdentity <- function(x, discard = 0, quiet = FALSE, plot = FALSE){
  
  tabname <- paste("acc", x@locus@sql, sep = "_")
  tabname <- gsub("__", "_", tabname)
  
  conn <- dbconnect(x)
  y <- paste("SELECT identity", "FROM", tabname)
  y <- dbGetQuery(conn, y)
  
  z <- paste("SELECT taxon, count(taxon), max(identity) AS maxident", 
             "FROM", tabname,
             "GROUP BY taxon", "ORDER BY maxident")
  z <- dbGetQuery(conn, z)
  dbDisconnect(conn)
  threshold <- z$maxident[discard + 1]
  
  ## screen output:
  ## --------------
  if ( !quiet ){
    discard <- ifelse(discard > 0, paste(" but", discard), "")
    out <- paste("threshold to include all", discard, 
                 " species at least once", sep = "")
    out <- c("locus", "minimum identity", "maximum identity", 
             "current threshold", out)
    out <- format(out, justify = "right")
    out[-1] <- paste("\n", out[-1], sep = "")
    val <- format(c(range(y$identity, na.rm = TRUE), 
                    x@locus@min.identity, threshold))
    cat(paste(out, c(x@locus@sql, val), sep = " : "), "\n")
  }
 
  if ( plot ){
    hist(y$identity, xlab = "Identity", main = x@locus@sql)
  } 
    
  list(current = x@locus@min.identity,
       suggested = threshold,
       ranking = z)
}