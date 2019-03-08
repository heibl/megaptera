## This code is part of the megaptera package 
##  © C. Heibl 2017 (last update 2019-03-05)

#' @title Utilities for NCBI Taxdump
#' @description Get subset of a taxonomy table in parent-child format.
#' @param tax Either an object of class \code{\link{megapteraProj}} or a data
#'   frame containing a taxonomy table in parent-child format as returned by
#'   \code{\link{dbReadTaxonomy}}.
#' @param mrca A character string giving the most recent common ancestor (MRCA)
#'   of the subset.
#' @param species A vector of mode \code{"character"} giving a subset of species
#'   names.
#' @param root A character string choosing between \code{"mrca"} and
#'   \code{"tol"} (tree of life).
#' @param syn Logical, indicating if synonyms should be returned in addition to
#'   accepted names; only has an effect if \code{tax} is of class
#'   \code{\link{megapteraProj}}.
#' @seealso \code{\link{dbReadTaxonomy}},
#'   \code{\link{taxdumpChildren}},\code{\link{taxdumpLineage}},
#'   \code{\link{taxdumpAddNode}}, \code{\link{taxdump2phylo}}.
#' @export

taxdumpSubset <- function(tax, mrca, species, root = "tol", syn = FALSE){
  
  root <- match.arg(root, c("tol", "mrca"))
  
  if (inherits(tax, "megapteraProj")){
    
    #############################################
    ## 'tax' is object of class 'megapteraProj'
    #############################################
    
    tax <- dbReadTaxonomy(tax, subset = species, syn = syn)
    
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
    
  } else {
    
    #############################################
    ## 'tax' is data frame in parent-child format
    #############################################
    
    # c("parent_id", "id", "taxon", "rank", "status")
    
    ## Get subset if species are given
    ## -------------------------------
    if (!missing(species)){
      
      ## Determine which separator is used by 'tax' 
      ## and impose it on 'species'
      ## --------------------------
      underscore <- length(grep("_", tax$taxon[tax$rank == "species"]))
      space <- length(grep(" ", tax$taxon[tax$rank == "species"]))
      if (underscore > space){
        species <- gsub(" ", "_", species)
      } else {
        species <- gsub("_", " ", species)
      }
      
      ## Check if 'tax' contains all 'tips'
      ## ----------------------------------
      not <- setdiff(species, tax$taxon)
      if (length(not)){
        cat("\nWARNING: ", length(not), " elements of 'tips' not in 'tax'", 
            paste0("\n- ", not), sep = "")
        species <- setdiff(species, not)
      }
      
      ## Create subset
      ## -------------
      id <- all_ids <- tax[tax$taxon %in% species, "id"]
      while (length(id) > 1) {
        id <- unique(tax$parent_id[tax$id %in% id & tax$status == "scientific name"])
        all_ids <- c(all_ids, id)
      }
      # tax <- tax[tax$id %in% all_ids, c("parent_id", "id", "taxon", "rank")]
      tax <- tax[tax$id %in% all_ids, ]
    }
    
    ## Delete tree-of-life root "tail" if desired
    ## ------------------------------------------
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
    } ## end of IF (line 87)
    
    ## Get subset if MRCA is given
    ## ---------------------------
    if (!missing(mrca)){
      id <- all_ids <- tax[tax$taxon %in% mrca, "id"]
      repeat {
        id <- tax[tax$parent_id %in% id, "id"]
        if (!length(id)) break
        all_ids <- c(all_ids, id)
      }
      ## add tree-of-life root "tail"
      if (root == "tol") all_ids <- unique(c(all_ids, taxdumpLineage(tax, mrca)$id))
      tax <- tax[tax$id %in% all_ids, ]
    }
  }
  tax
}
