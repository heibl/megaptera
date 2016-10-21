## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-01-20)

splitAlignment <- function(root, gt, a){
  root.descendants <- gt$edge[gt$edge[, 1] == root, 2]
  
  engine <- function(node, gt, a){
    
    if ( node > Ntip(gt) ){
      tips <- descendants(gt, node, labels = TRUE)
    } else {
      tips <-  gt$tip.label[node]
    }
    deleteEmptyCells(a[tips, ], quiet = TRUE)
    
  }
  lapply(root.descendants, engine, gt = gt, a = a)
}