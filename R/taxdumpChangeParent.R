## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-05-15)

#' @title Utilities for NCBI Taxdump
#' @description Change a node's parent in the NCBI taxonomy.
#' @param x An object of class \code{\link{megapteraProj}} or a \code{data.frame}
#'   a parent-child-format.
#' @param taxon A character string giving the name of the taxon.
#' @param new.parent A character string giving the name of the new parent taxon.
#' @return A \code{data.frame} in parent-child-format, if \code{class(x) ==
#'   "data.frame"}, or \code{TRUE}, if \code{class(x) == "megapteraProj"}. In
#'   the latter case \code{taxdumpChangeParent} is called for its side effect on
#'   the taxonomy table of the database.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpChildren}},
#'   \code{\link{taxdumpDropTip}}, \code{\link{taxdumpHigherRank}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
#' @export

taxdumpChangeParent <- function(x, taxon, new.parent){
  
  if (inherits(x, "megapteraProj")){
    tax <- dbReadTaxonomy(x)
  } else {
    tax <- x
  }
  
  if (!taxon %in% tax$taxon) {
    stop("Taxon '", taxon, "' not present in taxonomy", sep = "")
  }
  
  ## Look up IDs
  ## -----------
  opid <- tax$parent_id[tax$taxon == taxon]
  npid <- tax$id[tax$taxon == new.parent]

  
  if (inherits(x, "megapteraProj")){
    conn <- dbconnect(x)
    SQL <- paste("UPDATE taxonomy",
                 "SET", wrapSQL(npid, "parent_id", "="), 
                 "WHERE", wrapSQL(taxon, "taxon", "="))
    SQL <- dbGetQuery(conn, SQL)
    dbDisconnect(conn)
    cat("Parent of '", taxon, "'changed in taxonony table from '", tax$taxon[tax$id == opid],
        "' [", opid, "] to '", new.parent, "' [", npid, "]", sep = "")
    invisible(TRUE)
  } else {
    tax$parent_id[tax$taxon == taxon] <- npid
    return(tax)
  } 
}
