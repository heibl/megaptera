## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2015-06-10)

addRanks <- function(z, ranks){
  nr <- z$name[z$rank == "no rank"]
  id <- match(ranks, z$rank, incomparables = "no rank")
  zz <- z[id, ]
  zz$rank <- ranks
  
  nr <- z$name[z$rank == "no rank"]
  ff <- function(i, z){
    ff <- z$name[(which(z$name == i) + 1):nrow(z)]
    1:(min(which(z$name %in% ff)) - 1)
  }
  nrr <- lapply(nr, ff, z = z)
  names(nrr) <- nr
  nrl <- sapply(nrr, length)
  nid <- sort(nrl, FALSE)
  zz$name[nid] <- names(nid)
  
  ## fill remaining with "-"
  zz$name[is.na(zz$name)] <- "-"
  zz
}
