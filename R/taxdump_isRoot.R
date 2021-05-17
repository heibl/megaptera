## This code is part of the megaptera package
## © C. Heibl 2021 (last update 2021-05-17)

#' @title Utilities for NCBI Taxdump
#' @description Check if node is root node.
#' @param x A data frame containing a parent-child table as returned by
#'   \code{\link{dbReadTaxonomy}}.
#' @param taxon A character string giving the name of the node to be
#'   tested.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpAddNode}},
#'   \code{\link{taxdumpChildren}}, \code{\link{taxdumpDropTip}},
#'   \code{\link{taxdumpHigherRank}}, \code{\link{taxdumpMRCA}},
#'   \code{\link{taxdumpSubset}}, \code{\link{taxdump2phylo}}
#' @export

taxdump_isRoot <- function(x, taxon){
  
  x$taxon[x$id == x$parent_id[x$taxon == taxon]] == "root"
}