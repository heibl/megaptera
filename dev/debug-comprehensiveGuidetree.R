library(viper)

load("debugComprehensiveguidetree.rda")

tt <- subtrees[[3]]
bigtree(tt)
bigtree(fixNodes(tt), file = "nodesFixed.pdf")

id <- which(tt$edge[, 2] == 4811)
ttt <- tt
ttt$edge <- ttt$edge[-id, ]
ttt$Nnode <- ttt$Nnode - 1
bigtree(ttt, file = "edgeRemoved.pdf")


tre <- taxdumpDaughters(tax, "Adephaga", tip.rank = tip.rank)
edge.matrix <- as.matrix(tre[, c("parent_id", "id")])
rownames(edge.matrix) <- colnames(edge.matrix) <- NULL
## Sanity check for edge matrix

edgematrixSanity <- function(edge.matrix){
  
  sanity <- vector(length = 4)
  
  ## 1. Every node has only one parent node
  ##    i.e. there must not be duplicate node numbers in column 2
  sanity[1] <- !any(duplicated(edge.matrix[, 2]))
  
  ## 2. There can be only one root node
  ##    i.e. all but one nde numbers in colum 1 must also be 
  ##    present in column 2
  sanity[2] <- length(which(!edge.matrix[, 1] %in% edge.matrix[, 2])) == 1
  
  ## 3. The number of edges must the number of nodes - 1
  sanity[3] <- nrow(edge.matrix)  == length(unique(as.vector(edge.matrix))) - 1
  
  ## 4. There may not be loops, i.e. a node must not refer to itself
  sanity[4] <- !any(edge.matrix[, 1] == edge.matrix[, 2])
  all(sanity)
}

##
id <- which(!tre[, 2] %in% tre[, 1])
length(id)
test <- tre[id, ]


node <- 48622
is.terminal <- function(edge.matrix, node){
  ~
}
