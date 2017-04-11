## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-11-15)

#' @title Transitivity Merge
#' @description Merges two multiple sequence alignments using the transitivity 
#'   criterion.
#' @export

transitivityMerge <- function(obj, id = c(1, 2)){
  
  obj1 <- obj[[id[1]]]; obj2 <- obj[[id[2]]]
  
  ## shared species
  shared <- intersect(names(obj1), names(obj2))
  if ( length(shared) == 0 ) stop("MSAs do not overlap")
  
  ## identify all-gap positions
  ## --------------------------
  al1 <- max(sapply(obj1, max))
  some.bases <- unique(unlist(obj1[shared]))
  all.gap1 <- which(!(1:al1 %in% some.bases))
  
  al2 <- max(sapply(obj2, max))
  some.bases <- unique(unlist(obj2[shared]))
  all.gap2 <- which(!(1:al2 %in% some.bases))
  
  ##  step indices using all-gap information
  stepIndex <- function(index, all.gap.pos){
    for ( i in all.gap.pos ){
      index[index >= i ] <- index[index >= i ] + 1
    }
    index
  }
  obj1 <- lapply(obj1, stepIndex, all.gap.pos = all.gap2)
  obj2 <- lapply(obj2, stepIndex, all.gap.pos = all.gap1)
  
  ## remove one set of shared species
  obj1 <- obj1[!names(obj1) %in% shared]
  
  ## join bot sets and return
  c(obj1, obj2)
}