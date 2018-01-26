## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2018-01-15)

#' @title Utilities for NCBI Taxdump
#' @description Convert a taxonomy table in parent-child format into an object 
#'   of class \code{\link{phylo}}.
#' @param tax A data frame representing a taxonomy in parent-child format as 
#'   returned by \code{\link{dbReadTaxonomy}}.
#' @param tip.rank A character string giving the name of a rank. This rank will be
#'   treated as tip rank, i.e. all taxa of lower rank will be dicarded.
#' @return An object of class \code{\link{phylo}}.
#' @seealso \code{\link{dbReadTaxonomy}}, 
#'   \code{\link{taxdumpChildren}},\code{\link{taxdumpLineage}}, 
#'   \code{\link{taxdumpSubset}}, \code{\link{taxdumpAddNode}}.
#' @importFrom ape collapse.singles
#' @export

taxdump2phylo <- function(tax, tip.rank){
  
  tip.rank <- match.arg(tip.rank, c("species", "genus", "family"))
  
  ## Do some checks
  ## --------------
  tax <- unique(tax)
  if (is.character(tax$id)) tax$id <- as.numeric(tax$id)
  if (is.character(tax$parent_id)) tax$parent_id <- as.numeric(tax$parent_id)
  
  ## Truncate taxonomy to tip rank; this is done in two steps:
  ## ---------------------------------------------------------
  
  ## Step 1: Remove all nodes below tip.rank
  ## ---------------------------------------
  id <- tax[tax$rank == tip.rank, "id"]
  tdDescendants <- function(tax, id){
    all_ids <- vector()
    gain <- length(id)
    while (gain > 0){
      id <- tax[tax$parent_id %in% id, "id"]
      all_ids <- c(all_ids, id)
      gain <- length(id)
    }
    all_ids
  }
  id <- lapply(id, tdDescendants, tax = tax)
  id <- unlist(id)
  tax <- tax[!tax$id %in% id, ]
  
  ## Step 2
  ## There can be lineages with tip.rank missing,
  ## e.g. subgenus: Neocicindela, no genus, tribe: Cicindelini,
  ## These lineages will be dropped entirely.
  ## ----------------------------------------
  tn <- taxdump_isTerminal(tax)
  id <- tax$id[tn & tax$rank != tip.rank]
  if (length(id)){
    warning(length(id)," terminal taxa without a taxon of rank '", tip.rank, 
            "' in their lineage were removed:",
            paste("\n-", tax$taxon[tax$id %in% id]))
    tax <- taxdumpDropTip(tax, id)
  }
  
  ## BEGIN ASSEMBLING THE PHYLO OBJECT:
  ## ----------------------------------
  phy <- list()
  
  ## parent-child table is edge matrix
  ## ---------------------------------
  phy$edge <- as.matrix(tax[, c("parent_id", "id")])
  rownames(phy$edge) <- colnames(phy$edge) <- NULL
  storage.mode(phy$edge) <- "integer"
  
  ## SET TIP LABELS, RENUMBER TIPS
  ## -----------------------------
  tip.id <- tax$rank == tip.rank
  ntip <- length(which(tip.id)) 
  if (ntip == 1) return(tax$taxon[tip.id]) ## break here if only 1 leaf
  ## Make sure node numbers won't be overwritten
  phy$edge <- phy$edge + max(phy$edge) 
  phy$tip.label <- tax[tip.id, "taxon"]
  phy$tip.label <- gsub(" ", "_", phy$tip.label)
  phy$edge[tip.id, 2] <- seq_along(phy$tip.label)
  
  
  ## SET ROOT NODE
  ## -------------
  ## Identity root
  root.id <- setdiff(phy$edge[, 1], phy$edge[, 2])
  phy$edge <- phy$edge[phy$edge[, 2] != root.id, ] ## drop NCBI "root"
  ## If present, drop singleton root edge(s)
  while (length(phy$edge[phy$edge[, 1] == root.id, 2]) == 1){
    phy$edge <- phy$edge[phy$edge[, 1] != root.id, ] ## drop root edge
    root.id <- setdiff(phy$edge[, 1], phy$edge[, 2])
  }
  ## Renumber root node
  current.node <- length(phy$tip.label) + 1 ## will be reused below!
  phy$edge[phy$edge[, 1] == root.id, 1] <- current.node
  
  
  ## RENUMBER INTERNAL NODES
  ## -----------------------
  root_children <- phy$edge[phy$edge[, 1] == current.node, 2]
  alt_nodes <- root_children[-1]
  n <- root_children[1]
  while (length(alt_nodes) | length(n)){
    next_node <- phy$edge[phy$edge[, 1] == n, 2]
    if (length(next_node) == 0){
      ## Lineages ending before species rank:
      ## eliminate them
      delete_id <- phy$edge[, 1] == n
      phy$edge <- phy$edge[!delete_id, ]
      n <- tail(alt_nodes, 1)
      alt_nodes <- head(alt_nodes, -1)
      next
    }
    if (length(next_node) == 1){
      ## Delete singleton edges (e.g. monotypic genus)
      delete_id <- phy$edge[, 1] == n
      change_id <- phy$edge[, 2] == n
      phy$edge[change_id, 2] <- next_node
      phy$edge <- phy$edge[!delete_id, ]
      n <- next_node
      next
    }
    ## normal case: nodes of degree > 2
    id <- which(phy$edge == n, arr.ind = TRUE)
    current.node <- current.node + 1
    phy$edge[id] <- current.node
    next_node <- setdiff(next_node, 1:ntip)
    if (length(next_node)){
      n <- next_node[1]
      alt_nodes <- c(alt_nodes, next_node[-1])
    } else {
      n <- tail(alt_nodes, 1)
      alt_nodes <- head(alt_nodes, -1)
    }
  }

  
  
  ## REORDER EDGE MATRIX  
  ## -------------------
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
  ## ape::collapse.singles will not work when reorder does not
  ## work, which is often the case with taxonomies
  # phy <- collapse.singles(phy)
  phy$Nnode <- length(unique(phy$edge[, 1]))
  # phy <- collapseSingles(phy)
  phy
}
