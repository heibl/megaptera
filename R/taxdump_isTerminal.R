## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2018-01-15)

#' @title Utilities for NCBI Taxdump
#' @description Check if nodes are terminal, e.g. they are not parent nodes.
#' @param tax A data frame containing a parent-child table as returned by
#'   \code{\link{dbReadTaxonomy}}.
#' @param id A vector of mode \code{"numeric"} giving the ID of the nodes to be
#'   tested; can be missing an then all nodes are tested.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpAddNode}},
#'   \code{\link{taxdumpChildren}}, \code{\link{taxdumpDropTip}},
#'   \code{\link{taxdumpHigherRank}}, \code{\link{taxdumpMRCA}},
#'   \code{\link{taxdumpSubset}}, \code{\link{taxdump2phylo}}
#' @export

taxdump_isTerminal <- function(tax, id){
  
  ## Checks
  ## ------
  if (!is.data.frame(tax)){
    stop("'tax' is not a data frame")
  }
  if (!all(names(tax) %in% c("parent_id", "id", "taxon", "rank"))){
    stop("'tax' is not a valid parent-child table")
  }
  
  if (missing(id)) id <- tax$id
  !sapply(id, "%in%", table = tax[, "parent_id"])
}
