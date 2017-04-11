## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-02-20)

#' @export

splitAlignment <- function(root, gt, a){
  
  ## determine clades to be split
  root.descendants <- gt$edge[gt$edge[, 1] == root, 2]
  
  ## do the splitting
  engine <- function(node, gt, a){
    
    if ( node > Ntip(gt) ){
      tips <- descendants(gt, node, labels = TRUE)
    } else {
      tips <-  gt$tip.label[node]
    }
    deleteEmptyCells(a[tips, ], quiet = TRUE)
    
  }
  a <- lapply(root.descendants, engine, gt = gt, a = a)
  
  ## splitting can cause non-overlapping subalignments
  a <- lapply(a, splitNonoverlapping)
  
  ## unlist to first level
  unlistFirstLevel(a)
}