## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-02-22)

#' @export

taxdumpHigherRank <- function(x, taxon, rank){
  
  if (inherits(x, "megapteraProj")){
    x <- dbReadTaxonomy(x, subset = taxon)
  }
  x <- lapply(taxon, taxdumpLineage, x = x)
  sapply(x, function(z) z$taxon[z$rank == rank])
}

