## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-04-07)

#' @export

taxdumpMRCA <- function(x, species, tip.rank){
  
  ## Note: the tip.rank argument is necessary for cases
  ## when tip.rank = "genus", but the actual target is
  ## species level like in megaptera:stepD()
  
  if (!inherits(x, "megapteraProj")){
    stop("'x' is not of class 'megapteraProj'")
  }
  species <- gsub("_", " ", species)
  
  obj <- dbReadTaxonomy(x, tip.rank = tip.rank, subset = species, root = "mrca")
  obj[!obj$parent_id %in% obj$id, "taxon"]
}
