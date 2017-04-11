## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-01-12)

#' @export

splitGiTaxon <- function(x, enforce.binomial = FALSE, sep = " "){
  
  fun <- function(x){
    x <- unlist(strsplit(x, "_"))
    c(gi = tail(x, 1), 
      taxon = paste(head(x, -1), collapse = sep))
  }
  x <- lapply(x, fun)
  x <- do.call(rbind, x)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  
  if ( enforce.binomial ){
    x[, 2] <- strip.infraspec(x[, 2])
  }
  x
}