taxdumpDropTip <- function(x, id){
  
  ## Make sure that id is a terminal node
  ## ------------------------------------
  if (any(x$parent_id == id)) stop("'id' is not a terminal node")
  
  ## Move up the lineage until there is a sister node.
  ## All nodes from 'id' to the node that has a sister
  ## will be removed.
  ## ------------------------------------------------
  repeat {
    pid <- x[x$id == id[1], "parent_id"]
    n_children <- nrow(x[x$parent_id == pid, ])
    if (n_children > 1) break
    id <- c(pid, id)
  }
  
  x[!x$id %in% id, ]
}