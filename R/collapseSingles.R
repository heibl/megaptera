## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-11-20)

#' @title Collapse Singleton Edges
#' @description Collapse singleton edges in a phylogenetic tree. This function
#'   is an alternative for collapse.singleton from the \strong{ape} package that
#'   will also work when \code{ape::reorder(phy)} throws an error. (And this
#'   tends to happen often when taxonomic relationships are converted to trees.)
#' @param phy An object of class \code{"\link{phylo}"}.
#' @return An object of class \code{"\link{phylo}"}.
#' @seealso \code{\link{taxdump2phylo}}, \code{\link{comprehensiveGuidetree}}
#' @export

collapseSingles <- function(phy){
  
  ## Determine root node
  ## -------------------
  root <- setdiff(phy$edge[, 1], phy$edge[, 2])
  
  int_nodes <- setdiff(phy$edge[, 1], root)
  
  ## Helper function 1: Count number of descendents
  ## of a node
  ## ----------------------------------------------
  nb_descendents <- function(phy, node){
    length(phy$edge[phy$edge[, 1] == node, 2])
  }
  
  ## Helper function 2: Check if a node is a singleton,
  ## i.e. its descedents also have 0 (tips) or 1 descendent node
  ## -----------------------------------------------------------
  is_single <- function(phy, node){
    d <- phy$edge[phy$edge[, 1] == node, 2]
    d <- sapply(d, nb_descendents, phy = phy)
    all(d <= 1)
  }
  
  collapse <- sapply(int_nodes, is_single, phy = phy)
  collapse <- int_nodes[collapse]
  
  for (i in collapse){
    id <- phy$edge[, 2] == i
    subtending <- phy$edge[id, 1]
    phy$edge <- phy$edge[!id, ]
    phy$edge[phy$edge[, 1] == i, 1] <- subtending
  }
  
  ## Add number of internal nodes
  ## ----------------------------
  phy$Nnode <- max(phy$edge) - length(phy$tip.label)
  
  phy
  
  
}