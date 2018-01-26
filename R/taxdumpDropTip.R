## This code is part of the megaptera package 
##  © C. Heibl 2017 (last update 2018-01-15)

#' @title Utilities for NCBI Taxdump
#' @description Get subset of a taxonomy table in parent-child format.
#' @param tax Either an object of class \code{\link{megapteraProj}} or a data
#'   frame containing a taxonomy table in parent-child format as returned by
#'   \code{\link{dbReadTaxonomy}}.
#' @param id a vector of mode \code{"numeric"}, giving one or more IDs of a
#'   terminal node (tip).
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpAddNode}},
#'   \code{\link{taxdumpChildren}}, \code{\link{taxdumpHigherRank}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
#' @export


taxdumpDropTip <- function(tax, id){
  
  ## Make sure that id is a terminal node
  ## ------------------------------------
  isTerminal <- taxdump_isTerminal(tax, id)
  if (!all(isTerminal)) {
    stop("these nodes are not terminal: ", 
         paste(id[!isTerminal], collapse = ", "))
  }
  
  ## We have to loop over id because apply-like
  ## procedures would not handle sister lineages
  ## correctly.
  ## ----------
  for (i in id){
    this_id <- i
    ## Move up the lineage until there is a sister node.
    ## All nodes from 'id' to the node that has a sister
    ## will be removed.
    ## ------------------------------------------------
    repeat {
      pid <- tax[tax$id == this_id[1], "parent_id"]
      n_children <- nrow(tax[tax$parent_id == pid, ])
      if (n_children > 1) break
      this_id <- c(pid, this_id)
    }
    tax <- tax[!tax$id %in% this_id, ]
  }
  tax
}
