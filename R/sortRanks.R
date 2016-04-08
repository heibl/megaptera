## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-11-09)

sortRanks <- function(x){
  
  ranks <- c("superkingdom", "kingdom", "subkingdom",
             "phylum", "subphylum",
             "superclass","class", "subclass", "infraclass",
             "superorder", "order", "suborder", "infraorder", "parvorder",
             "superfamily", "family", "subfamily",
             "tribe", "subtribe",
             "genus", "subgenus",
             "species group", "species subgroup",
             "species")
  
  ## check for ranks not contained in 'ranks'
  ## ----------------------------------------
  unknown.ranks <- !unique(unlist(x)) %in% c("no rank", ranks)
  if ( any(unknown.ranks) ){
    stop("unknown rank: ", 
         unique(unlist(x))[unknown.ranks])
  }
  
  ## delete ranks from 'ranks', which are not
  ## contained in the actual data set
  ## --------------------------------
  absent.ranks <-  !ranks %in% unique(unlist(x))
  if ( any(absent.ranks) ){
    ranks <- ranks[!absent.ranks]
  }
  
  ## create unified rank set
  ## -----------------------
  r <- paste("^", ranks, "$", sep = "")
  unified <- lapply(x, grep, pattern = r[1])
  unified <- unique(unified)
  unified <- c(rep("no rank", max(unlist(unified)) - 1), 
               ranks[1])
  ## loop over each pair of neibouring ranks
  for ( i in 1:(length(ranks) - 1) ){
    gap <- paste(r[i:(i + 1)], collapse = "|")
    gap <- lapply(x, grep, pattern = gap)
    gap <- unique(gap)
    gap <- gap[sapply(gap, length) > 1]
    if ( length(gap) == 0 ){
      ## CASE 1: no element of 'x' has both ranks
      unified <- c(unified, ranks[i + 1])
    } else {
      ## CASE 2: at least 1 element of 'x' has both ranks
      gap <- sapply(gap, diff)
      unified <- c(unified, 
                   rep("no rank", max(unlist(gap)) - 1), 
                   ranks[i + 1])
    }
  }
  unified
}
