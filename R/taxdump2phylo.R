## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-06-12)

#' @title Utilities for NCBI Taxdump
#' @description Convert a taxonomy table in parent-child format into an object 
#'   of class \code{\link{phylo}}.
#' @param x A data frame representing a taxonomy in parent-child format as 
#'   returned by \code{\link{dbReadTaxonomy}}.
#' @param tip.rank A character string giving the name a rank. This rank will be
#'   treated as tip rank, i.e. all taxa of lower rank will be dicarded.
#' @return An object of class \code{\link{phylo}}.
#' @seealso \code{\link{dbReadTaxonomy}}, 
#'   \code{\link{taxdumpDaughters}},\code{\link{taxdumpLineage}}, 
#'   \code{\link{taxdumpSubset}}, \code{\link{taxdumpAddNode}}.
#' @export

taxdump2phylo <- function(x, tip.rank){
  
  tip.rank <- match.arg(tip.rank, c("species", "genus"))
  
  x <- unique(x)
  
  phy <- list()
  
  ## parent-child table is edge matrix
  ## ---------------------------------
  phy$edge <- as.matrix(x[, c("parent_id", "id")])
  rownames(phy$edge) <- colnames(phy$edge) <- NULL
  
  ## set tips
  ## --------
  tip.id <- x$rank == tip.rank
  ntip <- length(which(tip.id)) 
  if (ntip == 1) return(x$taxon[tip.id]) ## break here if only 1 leaf
  phy$edge <- phy$edge + ntip
  phy$tip.label <- x[tip.id, "taxon"]
  phy$tip.label <- gsub(" ", "_", phy$tip.label)
  phy$edge[tip.id, 2] <- seq_along(phy$tip.label)
  
  ## set root node
  ## -------------
  root.id <- setdiff(phy$edge[, 1], phy$edge[, 2])
  phy$edge <- phy$edge[phy$edge[, 2] != root.id, ]
  root.id <- which(phy$edge[, 1] == root.id)
  current.node <- length(phy$tip.label) + 1
  phy$edge[root.id, 1] <- current.node

  ## set internal nodes
  ## ------------------
  while (any(phy$edge[, 2] > current.node)){
  # while (current.node < 6238){
    n <- phy$edge[phy$edge[, 2] > current.node, 2][1]
    id <- which(phy$edge == n, arr.ind = TRUE)
    if (nrow(id) == 1){
      ## lineages ending before species rank:
      ## eliminate them
      phy$edge <- phy$edge[-id[, "row"], ]
    } else {
      if (nrow(id) == 2){
        ## singleton node (e.g. monotypic genus)
        id1 <- phy$edge[, 1] == n
        id2 <- phy$edge[, 2] == n
        phy$edge[id2, 2] <- phy$edge[id1, 2]
        phy$edge <- phy$edge[!id1, ]
      } else {
        ## normal case: nodes of degree > 2
        current.node <- current.node + 1
        phy$edge[id] <- current.node
      }
    }
  }
  # if (max(phy$edge) > 2 * ntip + 1) stop("setting internal nodes failed")
  ## reorder edges
  node <- ntip + 1
  nodes <- 1:nrow(phy$edge)
  for (i in nodes){
  # for (i in 1:8){
    # if (node <= ntip) next
    id <- which(phy$edge[, 1] == node)
    id <- c((1:i) - 1, ## head
            id, ## body
            setdiff(nodes, union((1:i) - 1, id)) ## tail
    )
    id <- id[id > 0]
    phy$edge <- phy$edge[id, ]
    node <- phy$edge[i, 2] ## next node
  }
  class(phy) <- "phylo"
  phy <- collapse.singles(phy)
  # test <- table(phy$edge[phy$edge > length(phy$tip.label)])
  # test <- as.numeric(names(test)[test == 1])
  # test.len <- length(test)


  # if (test.len > 0){
  #   phy$edge <- phy$edge[!(phy$edge[, 2] %in% test), ]
  #   test <- phy$edge[, 2][phy$edge[, 2] > min(test)]
  #   if (length(test) > 0){
  #     phy$edge[phy$edge %in% test] <- phy$edge[phy$edge %in% test] - (test.len - 1)
  #   }
  # }
  phy$Nnode <- max(phy$edge) - length(phy$tip.label)
  phy
}
