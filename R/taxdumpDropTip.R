## This code is part of the megaptera package 
##  © C. Heibl 2017 (last update 2017-11-16)

#' @title Utilities for NCBI Taxdump
#' @description Get subset of a taxonomy table in parent-child format.
#' @param tax Either an object of class \code{\link{megapteraProj}} or a data
#'   frame containing a taxonomy table in parent-child format as returned by
#'   \code{\link{dbReadTaxonomy}}.
#' @param id Numeric the ID of a terminal node (tip).
#' @seealso \code{\link{dbReadTaxonomy}},
#'   \code{\link{taxdumpDaughters}},\code{\link{taxdumpLineage}},
#'   \code{\link{taxdumpAddNode}}, \code{\link{taxdump2phylo}}.
#' @export


taxdumpDropTip <- function(tax, id){
  
  ## Make sure that id is a terminal node
  ## ------------------------------------
  if (any(tax$parent_id == id)) stop("'id' is not a terminal node")
  
  ## Move up the lineage until there is a sister node.
  ## All nodes from 'id' to the node that has a sister
  ## will be removed.
  ## ------------------------------------------------
  repeat {
    pid <- tax[tax$id == id[1], "parent_id"]
    n_children <- nrow(tax[tax$parent_id == pid, ])
    if (n_children > 1) break
    id <- c(pid, id)
  }
  
  tax[!tax$id %in% id, ]
}