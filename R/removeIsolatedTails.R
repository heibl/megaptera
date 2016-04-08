## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-02-25)

removeIsolatedTails <- function(x){
  
  core.function <- function(seqs){
    gaps <- seqs %in% as.raw(4)
    ## if no gaps are present, just return whole sequence
    if ( !any(gaps) ){
      return(seqs)
    }
    id <- vector(length = length(gaps) -1)
    for ( i in 1:(length(gaps) -1) )
      id[i] <- gaps[i] != gaps[i + 1]
    id <- which(id)
    id <- cbind(c(0, head(id, -1)) + 1, id)
    id <- cbind(id, id[, 2] - id[, 1] + 1)
    g <- if ( gaps[1] ) c(FALSE, TRUE) else c(TRUE, FALSE)
    g <- rep_len(g, nrow(id))
    id <- cbind(id, g)
    colnames(id) <- c("from", "to", "n", "gap")
    
    ## only one gap at the 3' end:
    if ( nrow(id) == 2 & id[2, "gap"] == 1 ){
      return(seqs)
    }
    
    check <- which(id[, "gap"] == 1)
    remove.until <- 0
    for (i in check ){
      if ( id[i, "n"] > id[i + 1, "n"] ){
        break
      } else {
        remove.until <- id[i + 1, "to"]
      }
    }
    if ( remove.until > 0 ) 
      seqs[1:remove.until] <- as.DNAbin(rep("-", remove.until))
    seqs
  }
  x <- apply(x, 1, core.function)
  x <- t(x)
  class(x) <- "DNAbin"
  x
}

