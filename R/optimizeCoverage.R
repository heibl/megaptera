## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2015-12-09)

optimizeCoverage <- function(x, discard = 0, quiet = FALSE, plot = FALSE){
  
  tabname <- paste("acc", x@locus@sql, sep = "_")
  tabname <- gsub("__", "_", tabname)
  
  conn <- dbconnect(x)
  y <- paste("SELECT coverage", "FROM", tabname)
  y <- dbGetQuery(conn, y)
  
  z <- paste("SELECT taxon, count(taxon), max(coverage) AS maxcov", 
             "FROM", tabname,
             "GROUP BY taxon", "ORDER BY maxcov")
  z <- dbGetQuery(conn, z)
  dbDisconnect(conn)
  threshold <- z$maxcov[discard + 1]
  
  ## screen output:
  ## --------------
  if ( !quiet ){  
    discard <- ifelse(discard > 0, paste(" but", discard), "")
    out <- paste("threshold to include all", discard, 
                 " species at least once", sep = "")
    out <- c("locus", "minimum coverage", "maximum coverage", 
             "current threshold", out)
    out <- format(out, justify = "right")
    out[-1] <- paste("\n", out[-1], sep = "")
    val <- format(c(range(y$coverage, na.rm = TRUE), 
                    x@locus@min.coverage, threshold))
    cat(paste(out, c(x@locus@sql, val), sep = " : "), "\n")
  }
  
  if ( plot ){
    hist(y$coverage, xlab = "Coverage", main = x@locus@sql)
  }
  
  list(current = x@locus@min.coverage,
       suggested = threshold,
       ranking = z)
}