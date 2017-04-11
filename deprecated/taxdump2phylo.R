## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-01-10)

#' @export

taxdump2phylo <- function(conn, taxon){
  
  x <- taxdumpDaughters(conn, taxon)
  
  id <- x$id[x$explode]
  while (length(id) > 0){
    x$explode <- FALSE
    for ( i in id){
      cat("\n", x$taxon[x$id == i])
      xx <- taxdumpDaughters(conn, i)
      xx$explode[xx$rank == "species"] <- FALSE
      if (!is.null(xx)){
        pos <- which(x$id == i)
        if (nrow(x) != pos){
          x <- rbind(x[1:pos, ], xx, x[(pos + 1):nrow(x), ])
        } else {
          x <- rbind(x, xx)
        }
      }
    } # end of FOR-loop
    id <- x$id[x$explode]
  } # end of WHILE-loop
  x$pid <- as.numeric(x$pid)
  x$id <- as.numeric(x$id)
  
  phy <- list()
  phy$edge <- as.matrix(x[, 1:2])
  rownames(phy$edge) <- colnames(phy$edge) <- NULL
  ## set tips
  tip.id <- x$rank == "species"
  ntip <- length(which(tip.id)) ## break here if only 1 leaf
  if (ntip == 1) return(x$taxon[tip.id])
  phy$edge <- phy$edge + ntip
  phy$tip.label <- x[tip.id, "taxon"]
  phy$edge[tip.id, 2] <- seq_along(phy$tip.label)
  ## set root node
  root.id <- phy$edge[, 1] == phy$edge[1, 1]
  current.node <- length(phy$tip.label) + 1
  phy$edge[root.id, 1] <- current.node
 
  ## set internal nodes
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
    message(current.node)
  }
  class(phy) <- "phylo"
  phy <- collapse.singles(phy)
  test <- table(phy$edge[phy$edge > length(phy$tip.label)])
  test <- as.numeric(names(test)[test == 1])
  test.len <- length(test)


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