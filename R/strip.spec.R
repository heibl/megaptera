## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-03-25)

#' @export

strip.spec <- function(x){
  if ( is.factor(x) ) x <- levels(x)[x]
  sepchar <- ifelse(length(grep("_", x)) != 0, "_", " ")
  x <- strsplit(x, sepchar)
  sapply(x, function(x) x[1])
}