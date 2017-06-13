## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-06-13)

#' @title Higher Taxa on Edges
#' @description Create a vector of taxon names to label edges in a phylogenetic 
#'   tree with \code{\link{edgelabels}}.
#' @param phy An object of class \code{\link{phylo}}.
#' @param tax A data frame containing a taxonomy table in parent-child format.
#' @return A data frame with two columns: 
#' \describe{ 
#'   \item{taxon}{the name of the taxon} 
#'   \item{edge}{the number of the edge in \code{phy}} 
#' }
#' @seealso \code{\link{dbReadTaxonomy}}
#' @export
#' @importFrom ape Nnode Ntip
#' @importFrom ips descendants

higherTaxaOnEdges <- function(phy, tax){
  int_nodes <- Ntip(phy) + (2:Nnode(phy))
  obj <- lapply(int_nodes, descendants, phy = phy, labels = TRUE)
  obj <- sapply(obj, taxdumpMRCA, x = tax)
  id <- which(phy$edge[, 2] %in% int_nodes)
  data.frame(taxon = obj, 
             edge = id,
             stringsAsFactors = FALSE)
}