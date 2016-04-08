## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-11-20)

splitGiTaxon <- function(x, enforce.binomial = FALSE){
  
  fun <- function(x){
    x <- unlist(strsplit(x, "_"))
    c(gi = tail(x, 1), 
      taxon = paste(head(x, -1), collapse = "_"))
  }
  x <- lapply(x, fun)
  x <- do.call(rbind, x)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  
  if ( enforce.binomial ){
    x[, 2] <- strip.infraspec(x[, 2])
  }
  x
}