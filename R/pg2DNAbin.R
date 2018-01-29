## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2018-01-26)

#' @export

pg2DNAbin <- function(pg){
  
  pg <- pg[order(pg$taxon, pg$pos), ]
  taxon <- unique(pg$taxon)
  
  ## Are sequences aligned
  aligned <- unique(tapply(pg$pos, pg$taxon, max))
  aligned <- ifelse(length(aligned) == 1, TRUE, FALSE)
  
  if (aligned){
    out <- matrix(NA, nrow = length(taxon), ncol = max(pg$pos))
    rownames(out) <- taxon
    for (i in 1:nrow(pg)){
      out[pg$taxon[i], pg$pos[i]] <- pg$nuc[i]
    }
  } else {
    out <- tapply(pg$nuc, pg$taxon, as.vector)
    names(out) <- taxon
  }
  
  as.DNAbin(out)
}
