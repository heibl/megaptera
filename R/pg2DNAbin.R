## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2018-01-30)

#' @export

pg2DNAbin <- function(pg){
  
  obj <- strsplit(pg$sequence, split = "")
  names(obj) <- pg$taxon
  obj <- as.DNAbin(obj)
  if (length(unique(sapply(obj, length))) == 1){
    obj <- as.matrix(obj)
  }
  obj
}
