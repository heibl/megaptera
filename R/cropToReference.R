## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-02-25)

cropToReference <- function(x){
  
  ## first and last position of reference
  ref <- x[grep("REF", rownames(x)), ]
  ref <- ref != as.raw(4)
  ref <- apply(ref, 1, function(z) range(which(z)))
  ref <- c(min(ref[1, ]), max(ref[, 2]))
  
  ## crop
  x[, ref[1]:ref[2]]
}