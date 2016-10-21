## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-01-20)

checkBlocks <- function(x, plot = TRUE, col){
  
  msa.tab <- dbTableNames(x, x@taxon@tip.rank)

  ## update spec_* table
  ## -------------------
  conn <- dbconnect(x)
  SQL <- paste("SELECT status, count(status)", 
               "FROM",  msa.tab,
               "GROUP BY status")
  b <- lapply(SQL, dbGetQuery, conn = conn)
  dbDisconnect(conn)
  names(b) <- msa.tab
  
  b <- b[sapply(b, nrow) > 0]
  b <- lapply(b, function(z) z$count)
  b <- lapply(b, sort, decreasing = TRUE)
  b <- b[order(sapply(b, head, n = 1))]
  
  mb <- max(sapply(b, length))
  foo <- function(z, mb){
    if ( length(z)  < mb )
      c(z, rep(0, mb - length(z)))
    else
      z
  }
  bb <- lapply(b, foo, mb = mb)
  names(bb) <- gsub(paste(x@taxon@tip.rank, "_", sep = ""), "", names(bb))
  names(bb) <- gsub("_", ".", names(bb))
  if ( plot ) {
    if ( missing(col) ) col <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
    barplot(t(do.call(rbind, bb)), 
                      horiz = TRUE, las = 1, col = col)
  }
  b
}