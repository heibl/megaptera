## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2018-01-31)

#' @title Summarize and Plot Aligment Blocks
#' @description Summarize and visualize alignment blocks with a barplot.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param plot Logical, for suppressing the barplot.
#' @param col A vector of mode \code{"character"} giving two or more colors two
#'   distinguish between alignment blocks in the barplot.
#' @param subset A vector of mode \code{"character"}, that can be used to choose
#'   a subset of the total taxa available.
#' @seealso \code{\link{checkSpecLocus}}
#' @importFrom DBI dbDisconnect dbGetQuery
#' @importFrom graphics barplot
#' @export

checkBlocks <- function(x, plot = TRUE, col, subset){
  
  if (!missing(subset)){
    subset <- gsub("_", " ", subset)
    subset <- paste0("^", subset, "$")
    subset <- paste(subset, collapse = "|")
    subset <- wrapSQL(subset)
    subset <- paste("AND", subset)
    
  } else {
    subset <- ""
  }
  
  ## How many blocks per locus?
  ## --------------------------
  conn <- dbconnect(x)
  loci <- "SELECT DISTINCT locus FROM species_sequence ORDER BY locus"
  loci <- dbGetQuery(conn, loci)$locus
  
  b <- paste("SELECT status, count(status)", 
             "FROM species_sequence", 
             "WHERE", wrapSQL(loci, "locus", "=", NULL),
             subset,
             "GROUP BY status")
  b <- lapply(b, dbGetQuery, conn = conn)
  dbDisconnect(conn)
  names(b) <- loci
  # b <- b[sapply(b, nrow) > 0]
  b <- lapply(b, function(z) z$count)
  b <- lapply(b, sort, decreasing = TRUE)
  b <- b[order(sapply(b, head, n = 1))]
  
  
  ## Optional: produce barplot
  ## -------------------------
  if (plot) {
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
    
    if (missing(col)) col <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
    barplot(t(do.call(rbind, bb)), 
            horiz = TRUE, las = 1, col = col)
  }
  
  ## Return list of blocks per locus
  ## -------------------------------
  b
}