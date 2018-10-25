## This code is part of the megaptera package 
##  © C. Heibl 2017 (last update 2018-06-26)

#' @title Utilities for NCBI Taxdump
#' @description Convert a species list into a taxonomy table in parent-child format.
#' @param species A vector of mode \code{"character"} giving subset.
#' @seealso \code{\link{dbReadTaxonomy}},
#'   \code{\link{taxdumpChildren}},\code{\link{taxdumpLineage}},
#'   \code{\link{taxdumpAddNode}}, \code{\link{taxdump2phylo}}.
#' @export

speclist2taxdump <- function(species){

  tax <- data.frame(id = 1, parent_id = 1, taxon = "root", ## 'root' required by taxdumpLineage! 
                    rank = "no rank", stringsAsFactors = FALSE)
  genera <- unique(strip.spec(species))
  genera <- data.frame(id = seq_along(genera) + 1, parent_id = 1, taxon = genera, rank = "genus", 
                       stringsAsFactors = FALSE)
  species <- data.frame(id = seq_along(species) + max(genera$id), parent_id = NA, 
                        taxon = species, rank = "species", stringsAsFactors = FALSE)
  tax <- rbind(tax, genera, species)
  for (i in which(is.na(tax$parent_id))){
    tax$parent_id[i] <- tax$id[tax$taxon == strip.spec(tax$taxon[i])]
  }
  tax
}