## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2018-12-20)

#' @title Utilities for NCBI Taxdump
#' @description Get the taxon name of a certain higher rank for a particular
#'   taxon, e.g. get the familiy of \emph{Megaptera}.
#' @param x Either an object of class \code{\link{megapteraProj}} or a data
#'   frame containing a taxonomy table in parent-child format as returned by
#'   \code{\link{dbReadTaxonomy}}.
#' @param taxon A vector of mode \code{"character"} giving the taxon name.
#' @param rank A vector of mode \code{"character"} giving the name of a higher
#'   (relative to \code{taxon}) rank.
#' @param it.max An integer giving maximum number of iteration in case
#'   \code{rank} is not present in the lineage of \code{taxon}
#' @return A character string giving a taxon name.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpAddNode}},
#'   \code{\link{taxdumpChildren}}, \code{\link{taxdumpDropTip}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
#' @export

taxdumpHigherRank <- function(x, taxon, rank, it.max = 30){
  
  if (inherits(x, "megapteraProj")){
    x <- dbReadTaxonomy(x, subset = taxon)
  }
  ## Algorithmus  1
  # x <- lapply(taxon, taxdumpLineage, tax = x, highest.rank = rank)
  # sapply(x, function(z) z$taxon[z$rank == rank]) # or
  # sapply(x, function(z) tail(z$taxon, 1)) # not faster!
  
  ## Algorithmus 2 is a bit faster
  core <- function(x, taxon, rank, it.max){
    z <- x[x$taxon == taxon & x$status == "scientific name", ]
    i <- 1
    repeat {
      i <- i + 1
      z <- x[x$id %in% x$parent_id[x$id == z$id] & x$status == "scientific name", ]
      if (z$rank == rank | i > it.max) break
    }
    ifelse(z$rank == rank, z$taxon, NA)
  }
  sapply(taxon, core, x = x, rank = rank, it.max = it.max)
}

