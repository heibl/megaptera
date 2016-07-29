## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-07-29)

findRoot <- function(tax){
  
  id <- apply(tax, 2, function(y) length(unique(y)))
  id <- cbind(cumsum(id), cumsum(rep(1, length(id))))
  id <- max(which(id[, 1] == id[, 2]))
  tax[1, id]
}