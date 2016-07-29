## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-07-27)

sortRanks <- function(x){
  
  ## these are the ordered build-in ranks
  ## ------------------------------------
  ranks <- c("superkingdom", "kingdom", "subkingdom",
             "phylum", "subphylum",
             "superclass","class", "subclass", "infraclass",
             "superorder", "order", "suborder", "infraorder", "parvorder",
             "superfamily", "family", "subfamily",
             "tribe", "subtribe",
             "genus", "subgenus",
             "species group", "species subgroup",
             "species")
  
  ## get the unique set of rank combinations
  ## and return immediately if set is of length 1
  ## --------------------------------------------
  rank.set <- unique(lapply(x, function(z) z$rank))
  if ( length(rank.set) == 1 ){
    cat("\n.. unique rank set ..")
    return(rank.set[[1]])
  }
  
  ## check for ranks not contained in 'ranks'
  ## ----------------------------------------
  unknown.ranks <- !unique(unlist(rank.set)) %in% c("no rank", ranks)
  if ( any(unknown.ranks) ){
    stop("unknown rank: ", 
         unique(unlist(rank.set))[unknown.ranks])
  }
  
  ## delete ranks from 'ranks', which are not
  ## contained in the actual data set
  ## --------------------------------
  absent.ranks <-  !ranks %in% unique(unlist(rank.set))
  if ( any(absent.ranks) ){
    ranks <- ranks[!absent.ranks]
  }
  
  ## create unified rank set
  ## -----------------------
  r <- paste("^", ranks, "$", sep = "")
  
  ## 1: check if there is a 'no rank' above highest rank
  ## ---------------------------------------------------
  unified <- lapply(rank.set, grep, pattern = r[1])
  unified <- unique(unified)
  unified <- c(rep("no rank", max(unlist(unified)) - 1), 
               ranks[1])
  
  
  ## loop over each pair of neighbouring ranks
  for ( i in 1:(length(ranks) - 1) ){
    
    ## Calculate number of gaps between elements
    ## i (A), i + 1 (B) and i + 2 (C);
    ## 0 means no gap and -1 means no element in 
    ## rank.set contains both ranks
    gap1 <- paste(r[i + (0:1)], collapse = "|")
    gap1 <- lapply(rank.set, grep, pattern = gap1)
    gap1 <- lapply(gap1, diff)
    gap1[sapply(gap1, length) == 0] <- 0
    
    gap2 <- paste(r[c(i, i + 2)], collapse = "|")
    gap2 <- lapply(rank.set, grep, pattern = gap2)
    gap2 <- lapply(gap2, diff)
    gap2[sapply(gap2, length) == 0] <- 0
    
    gap3 <- paste(r[i + (1:2)], collapse = "|")
    gap3 <- lapply(rank.set, grep, pattern = gap3)
    gap3 <- lapply(gap3, diff)
    gap3[sapply(gap3, length) == 0] <- 0
    
    gap <- data.frame(AB = unlist(gap1), 
                      AC = unlist(gap2), 
                      BC = unlist(gap3))
    gap <- unique(gap) - 1
    n <- apply(gap, 2, max)
    
    
    if ( n[1] == -1 ){
      ## CASE 1: no element of 'rank.set' has both ranks
      unified <- c(unified, ranks[i + 1])
    } else {
      ## CASE 2: at least 1 element of 'rank.set' has both ranks
      if ( n["BC"] == -1 &  n["AC"] > 0){
        unified <- c(unified, 
                     rep("no rank", n["AB"]), 
                     ranks[i + 1],
                     rep("no rank", n["AC"]))
      } else {
        unified <- c(unified, 
                     rep("no rank", n["AB"]), 
                     ranks[i + 1])
      }
    }
  }
  unified
}
