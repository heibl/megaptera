## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2021-03-15)

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
  
  ## Get taxonomy
  if (inherits(x, "megapteraProj")){
    x <- dbReadTaxonomy(x)
  }
  taxon <- gsub("_", " ", taxon)
  
  ## Make sure 'rank' is available in 'tax'
  rank <- match.arg(rank, unique(x$rank))
  
  ## Make sure 'taxon' is available in 'tax'
  if (!taxon %in% x$taxon){
    warning("there is no taxon '", taxon, "' in 'x'")
    return(NA)
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
      if (!nrow(z)) return(NA) ## footnote 1
      if (z$rank == rank | i > it.max) break
    }
    ifelse(z$rank == rank, z$taxon, NA)
  }
  sapply(taxon, core, x = x, rank = rank, it.max = it.max)
}

## footnote 1: This is for cases when the target rank is not present in the
## taxonomy and the 'tree-of-life tail' is shorter than it.max

