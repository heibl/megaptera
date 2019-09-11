## This code is part of the megaptera package 
##  © C. Heibl 2017 (last update 2019-05-03)

#' @title Utilities for NCBI Taxdump
#' @description Convert a list of species or genera into a taxonomy table in
#'   parent-child format.
#' @param species A vector of mode \code{"character"} giving a set of species or
#'   genus names.
#' @seealso \code{\link{dbReadTaxonomy}},
#'   \code{\link{taxdumpChildren}},\code{\link{taxdumpLineage}},
#'   \code{\link{taxdumpAddNode}}, \code{\link{taxdump2phylo}}.
#' @export

speclist2taxdump <- function(species){

  tax <- data.frame(id = 1, parent_id = 1, taxon = "root", ## 'root' required by taxdumpLineage! 
                    rank = "no rank", status = "scientific name",
                    stringsAsFactors = FALSE)
  genera <- strip.spec(species)
  tip.rank <- ifelse(all(genera == species), "genus", "species")
  
  genera <- unique(genera)
  genera <- data.frame(id = seq_along(genera) + 1, parent_id = 1, taxon = genera, 
                       rank = "genus", status = "scientific name",
                       stringsAsFactors = FALSE)
  species <- data.frame(id = seq_along(species) + max(genera$id), parent_id = NA, 
                        taxon = species, rank = "species", status = "scientific name",
                        stringsAsFactors = FALSE)
  if (tip.rank == "species") {
    tax <- rbind(tax, genera, species)
  } else {
    tax <- rbind(tax, genera)
  }
  for (i in which(is.na(tax$parent_id))){
    tax$parent_id[i] <- tax$id[tax$taxon == strip.spec(tax$taxon[i])]
  }
  tax
}