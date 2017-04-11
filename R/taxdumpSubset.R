## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-03-24)

#' @export

taxdumpSubset <- function(tax, mrca, species){
  
  if (inherits(tax, "megapteraProj")){
    tax <- dbReadTaxonomy(tax, subset = species)
  } else {
    ## get subset if tips are given
    if (!missing(species)){
      id <- all_ids <- tax[tax$taxon %in% species, "id"]
      while (length(id) > 1) {
        id <- unique(tax[tax$id %in% id, "parent_id"])
        all_ids <- c(all_ids, id)
      }
    }
    ## get subset if MRCA is given
    if (!missing(mrca)){
      id <- all_ids <- tax[tax$taxon %in% mrca, "id"]
      repeat {
        id <- tax[tax$parent_id %in% id, "id"]
        if (!length(id)) break
        all_ids <- c(all_ids, id)
      }
      ## add root "tail"
      all_ids <- unique(c(all_ids, taxdumpLineage(tax, mrca)$id))
    }
    tax <- tax[tax$id %in% all_ids, c("parent_id", "id", "taxon", "rank")]
  }
  tax
}