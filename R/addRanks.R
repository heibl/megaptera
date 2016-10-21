## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2015-07-27)

## Calathus vivesi: Calathus is included twice in the lineage
## as genus and as no rank. I included (2016-07-26) a patch that removes
## the no-rank-row in these cases.

addRanks <- function(z, ranks){
  
  ## tackle the 'Calathus' issue (see above)
  ## ---------------------------------------
  dups <- table(z$name)
  if ( any(dups > 1) ){
    id <- which(z$name == names(dups)[dups > 1])
    dr <- z$rank[id]
    cat("\nduplicate names in lineage:", z$name[z$rank == "species"])
    if ( all(dr %in% c("genus", "no rank")) ){
      z <- z[-id[dr == "no rank"], ]
    } else {
      # stop("Calathus-issue-type-error in addRanks")
    }
  }
  ## this is for debugging:
  # which(sapply(y, function(x) any(x$name %in% "Aphodius pedellus")))
  
  ## extend dataframe to complete lineage, i.e. include ranks 
  ## that are given for this lineage
  ## -------------------------------
  id <- match(ranks, z$rank, incomparables = "no rank")
  cl <- z[id, ]
  cl$rank <- ranks
  cl$name[is.na(cl$name) & cl$rank != "no rank"] <- "-"
  
  nr <- z$name[z$rank == "no rank"]
  ff <- function(i, z){
    id <- which(z$name == i)
    ff <- z$name[(id + 1):nrow(z)]
    1:(min(which(z$name %in% ff)) - 1)
  }
  nrr <- lapply(nr, ff, z = z)
  names(nrr) <- nr
  nrl <- sapply(nrr, length)
  nid <- sort(nrl, FALSE)
  cl$name[nid] <- names(nid)
  
  ## fill remaining with "-"
  cl$name[is.na(cl$name)] <- "-"
  cl
}
