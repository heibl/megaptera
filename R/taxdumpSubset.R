## This code is part of the megaptera package © C. Heibl 2017 (last update
## 2017-06-12)

#' @title Utilities for NCBI Taxdump
#' @description Get subset of a taxonomy table in parent-child format.
#' @param tax Either an object of class \code{\link{megapteraProj}} or a data
#'   frame containing a taxonomy table in parent-child format as returned by \code{\link{dbReadTaxonomy}}.
#' @param mrca A character string giving the most revent common ancestor (mrca) of the subset.
#' @param species A vector of mode \code{"character"} giving subset.
#' @param root A character string choosing between \code{"mrca"} and
#'   \code{"tol"} (tree of life).
#' @seealso \code{\link{dbReadTaxonomy}},
#'   \code{\link{taxdumpDaughters}},\code{\link{taxdumpLineage}},
#'   \code{\link{taxdumpAddNode}}, \code{\link{taxdump2phylo}}.
#' @export

taxdumpSubset <- function(tax, mrca, species, root = "tol"){
  
  root <- match.arg(root, c("tol", "mrca"))
  
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
      tax <- tax[tax$id %in% all_ids, c("parent_id", "id", "taxon", "rank")]
      ## delete tree-of-life root "tail"
      if (root == "mrca") {
        id <- tax$id[tax$taxon == "cellular organisms"]
        repeat {
          new_id <- tax[tax$parent_id == id[1], "id"]
          if (length(new_id) > 1) {
            id[1] <- 1 ## replace MRCA by root
            break
          }
          id <- c(new_id, id)
        }
        tax <- tax[!tax$id %in% id, ]
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
      ## add tree-of-life root "tail"
      if (root == "tol") all_ids <- unique(c(all_ids, taxdumpLineage(tax, mrca)$id))
      tax <- tax[tax$id %in% all_ids, c("parent_id", "id", "taxon", "rank")]
    }
  }
  tax
}