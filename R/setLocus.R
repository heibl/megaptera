## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-10-28)

#' @export

setLocus <- function(x, locus){
  #   stopifnot(inherits(x, "megapteraProj"))
  #   stopifnot(inherits(locus, "locus"))
  x@locus <- locus
  x
}