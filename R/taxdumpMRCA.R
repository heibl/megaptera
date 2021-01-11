## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2019-11-13)

#' @title Utilities for NCBI Taxdump
#' @description Get most-recent common ancestor (MRCA) of a group of species.
#' @param x Either an object of class \code{\link{megapteraProj}} or a data
#'   frame containing a taxonomy table in parent-child format as returned by
#'   \code{\link{dbReadTaxonomy}}.
#' @param species A vector of mode \code{"character"} giving the species names.
#'   Can be missing, in which case the root is returned.
#' @param tip.rank A character string giving the name a rank.
#' @return A character string giving the name a the MRCA.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpAddNode}},
#'   \code{\link{taxdumpChildren}}, \code{\link{taxdumpDropTip}},
#'   \code{\link{taxdumpHigherRank}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
#' @export

taxdumpMRCA <- function(x, species, tip.rank){
  
  ## Note: the tip.rank argument is necessary for cases
  ## when tip.rank = "genus", but the actual target is
  ## species level like in megaptera:stepD()
  
  ## Create subset of taxonomy ...
  ## -----------------------------
  if (!missing(species)){
    
    ## species must not use underscores!
    species <- gsub("_", " ", species)
    
    if (inherits(x, "megapteraProj")){
      x <- dbReadTaxonomy(x, tip.rank = tip.rank, subset = species, root = "mrca")
    } else {
      x <- taxdumpSubset(x, species = species, root = "mrca")
    }
    
    ## ... or use whole taxonomy
    ## -------------------------
  } else {
    if (inherits(x, "megapteraProj")){
      x <- dbReadTaxonomy(x, tip.rank = tip.rank, root = "mrca")
    }
  }
  
  ## Identify MRCA
  ## The loop is necessary because 'x' can be with root='tol'
  ## --------------------------------------------------------
  x <- x[x$id != 1, ]
  id <- x[!x$parent_id %in% x$id, "id"]
  repeat {
    temp <- x$id[x$parent_id == id]
    if (length(temp) == 1){
      id <- temp
    } else {
      break
    }
  }
  x[x$id == id, "taxon"]
}
