## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-06-12)

#' @title Utilities for NCBI Taxdump
#' @description Get most-recent common ancestor (MRCA) of a group of species.
#' @param x Either an object of class \code{\link{megapteraProj}} or a data
#'   frame containing a taxonomy table in parent-child format as returned by \code{\link{dbReadTaxonomy}}.
#' @param species A vector of mode \code{"character"} giving the species names.
#' @param tip.rank A character string giving the name a rank.
#' @return A character string giving the name a the MRCA.
#' @seealso \code{\link{dbReadTaxonomy}},
#'   \code{\link{taxdumpDaughters}},\code{\link{taxdumpLineage}},
#'   \code{\link{taxdumpAddNode}}, \code{\link{taxdump2phylo}}.
#' @export

taxdumpMRCA <- function(x, species, tip.rank){
  
  ## Note: the tip.rank argument is necessary for cases
  ## when tip.rank = "genus", but the actual target is
  ## species level like in megaptera:stepD()
  
  ## species must not use underscores!
  species <- gsub("_", " ", species)
  
  if (inherits(x, "megapteraProj")){
    x <- dbReadTaxonomy(x, tip.rank = tip.rank, subset = species, root = "mrca")
  } else {
    x <- taxdumpSubset(x, species = species, root = "mrca")
  }
  x[!x$parent_id %in% x$id, "taxon"]
}
