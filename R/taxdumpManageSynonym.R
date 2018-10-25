## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2018-10-25)

#' @title Utilities for NCBI Taxdump
#' @description Manage Synonyms
#' @param tax A data frame in parent-child format.
#' @param binomials A vector of mode \code{"character"} giving one or more taxon
#'   names. The first name is considered to be the accepted name, the.
#' @param keep.syn Logical, should synonyms present in \code{tax} be kept?
#' @param add.syn Logical, should synonyms not present in \code{tax} be added?
#' @param quiet Logical, if \code{FALSE}, progress is printed on screen for debugging.
#' @seealso \code{\link{taxdumpAddNode}}, \code{\link{taxdumpChildren}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}}
#' @export

taxdumpManageSynonym <- function(tax, binomials, keep.syn = TRUE, add.syn = FALSE, quiet = TRUE){
  
  ## Assume that first element of 'binomials' is
  ## the accepted name, all others are synonms
  
  if (!quiet){
    cat("Accepted name:", binomials[1], ".. ")
  }
  id <- which(binomials %in% tax$taxon)
  if (!length(id)){
    cat("WARNING: not present (incl. synonyms)\n")
    return(tax)
  }
  
  ## 1. Add accepted name (if it is missing)
  ## ---------------------------------------
  if (min(id) > 1){

    ## Insert missing accepted name
    ## ----------------------------
    parent <- strip.spec(binomials[1])
    if (!parent %in% tax$taxon){
      parent <- sapply(binomials[id], taxdumpHigherRank, x = tax, rank = "family")
      ## DIRTY, DIRTY!
      parent <- table(parent)
      parent <- names(parent)[which.max(parent)]
    } 
    tax <- taxdumpAddNode(taxon = binomials[1], x = tax, rank = "species", parent = parent)
    tax$id[tax$taxon == binomials[1]]
    cat("ADDED AS", parent, "\n")
  } else {
    tax$status[tax$taxon == binomials[1]] <- "scientific name"
    cat("OK\n")
  }
  
  ## 2. Set ID of synonyms to ID of accepted name
  ## --------------------------------------------
  tax$id[tax$taxon %in% binomials[id]] <- tax$id[tax$taxon == binomials[1]]
  
  ## 3. Manage synonyms
  ## ------------------
  if (keep.syn){
    
    ## Set status to "synonym"
    tax$status[tax$taxon %in% binomials[id[id > 1]]] <- "synonym"
    
    if (add.syn){
     syn_to_add <- binomials[-id]
     if (length(syn_to_add)){
       family <- sapply(binomials[1], taxdumpHigherRank, x = tax, rank = "family")
       for (i in syn_to_add){
         tax <- taxdumpAddNode(taxon = i, x = tax, rank = "species", parent = family)
         ## CAUTION: Will set status to "scientific name"
       }
     }
    }
  } else {
    
    ## Delete synonyms (if there are any and 'keep.syn' == FALSE)
    if (length(id[id > 1])){
      tax <- tax[!tax$taxon %in% binomials[id[id > 1]], ]
    }
  }
  tax
}