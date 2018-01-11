## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2018-01-10)

#' @title Utilities for NCBI Taxdump
#' @description Get the taxon name of a certain higher rank for a particular
#'   taxon, e.g. get the familiy of \emph{Megaptera}.
#' @param x Either an object of class \code{\link{megapteraProj}} or a data
#'   frame containing a taxonomy table in parent-child format as returned by
#'   \code{\link{dbReadTaxonomy}}.
#' @param taxon A vector of mode \code{"character"} giving the taxon name.
#' @param rank A vector of mode \code{"character"} giving the name of a higher (relative to
#'   \code{taxon}) rank.
#' @return A character string giving a taxon name.
#' @seealso \code{\link{dbReadTaxonomy}},
#'   \code{\link{taxdumpDaughters}},\code{\link{taxdumpLineage}},
#'   \code{\link{taxdumpAddNode}}, \code{\link{taxdump2phylo}}.
#' @export

taxdumpHigherRank <- function(x, taxon, rank){
  
  if (inherits(x, "megapteraProj")){
    x <- dbReadTaxonomy(x, subset = taxon)
  }
  x <- lapply(taxon, taxdumpLineage, tax = x)
  sapply(x, function(z) z$taxon[z$rank == rank])
}

