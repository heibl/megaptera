## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-02-20)

#' @export

splitNonoverlapping <- function(a){
  
  if (nrow(a) == 1) return(a)
  if (nrow(a) > 2){
    # stop("more than 2 sequences in splitNonoverlapping")
    return(a)
  }
  engine <- function(z){
    which(z %in% as.raw(c(136, 40, 72, 24)))
  }
  id <- apply(a, 1, engine)
  if (is.matrix(id)) id <- list(id[, 1], id[, 2])
  if (length(intersect(id[[1]], id[[2]])) == 0){
    a <- list(a[1, ], a[2, ])
    a <- lapply(a, deleteEmptyCells, quiet = TRUE)
  }
  a
}