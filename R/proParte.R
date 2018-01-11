## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-11-06)

#' @keywords internal
#' importFrom ips sister
#' export

proParte <- function(phy, id){
  
  s1 <- lapply(id, sister, phy = phy)
  s2 <- sapply(s1, function(x, id) all(x %in% id), id = id)
  pp <- as.list(id[!s2])
  
  if ( any(s2) ){
    
    identifyBreaks <- function(z){
      
      zz <- list()
      while (length(z)){
        y <- z == (seq_len(length(z)) + min(z) - 1)
        zz <- c(zz, list(z[y]))
        z <- z[!y]
      }
      zz
    }
    s3 <- union(id[s2], unlist(s1[s2]))
    s3 <- identifyBreaks(s3)
    pp <- c(pp, s3)
  }
  pp
}