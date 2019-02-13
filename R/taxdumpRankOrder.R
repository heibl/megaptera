## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2018-12-19)

#' @title Utilities for NCBI Taxdump
#' @description Return the order of ranks in a taxonomy from highest to lowest.
#' @param x Either an object of class \code{\link{megapteraProj}} or a data
#'   frame containing a taxonomy table in parent-child format as returned by
#'   \code{\link{dbReadTaxonomy}}.
#' @return A vector of mode \code{"character"} containing the ordered ranks.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpAddNode}},
#'   \code{\link{taxdumpChildren}}, \code{\link{taxdumpDropTip}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
#' @export

taxdumpRankOrder <- function(x){
  
  if (inherits(x, "megapteraProj")){
    x <- dbReadTaxonomy(x, subset = taxon)
  }
  
  ranks <- unique(x$rank)
  ranks <- ranks[ranks != "no rank"]
  
  parentRanks <- function(x, r){
    unique(x$rank[x$id %in% x$parent_id[x$rank == r]])
  }
  p <- lapply(ranks, parentRanks, x = x)
  names(p) <- ranks
  
  for (i in seq_along(p)){
    has_parent <- sapply(p, function(a, b) any(a %in% b), b = names(p)[i])
    p[has_parent] <- lapply(p[has_parent], union, y = p[[i]])
  }
  id <- sapply(p, length)
  names(sort(id))
}