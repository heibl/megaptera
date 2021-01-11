## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2020-02-20)

#' @export

rcString <- function(s, reverse = TRUE, complement = TRUE){
  
  custom_rev <- function(z){
    paste(rev(unlist(strsplit(z, ""))), collapse = "")
  }
  
  if (complement){
    s <- gsub("a", "G", s)
    s <- gsub("g", "A", s)
    s <- gsub("c", "T", s)
    s <- gsub("t", "C", s)
    s <- gsub("r", "Y", s) # R = A|G
    s <- gsub("y", "R", s) # Y = C|T
    s <- gsub("w", "S", s) # W = A|T
    s <- gsub("s", "W", s) # S = G|C
    s <- gsub("m", "K", s) # M = A|C
    s <- gsub("k", "M", s) # K = G|T
    s <- gsub("h", "B", s) # H = A|T|C != G
    s <- gsub("b", "H", s) # B = G|C|T != A
    s <- gsub("v", "D", s) # V = G|A|C != T
    s <- gsub("d", "V", s) # D = A|G|T != C
  }
  
  if (reverse){
    s <- sapply(s, custom_rev)
  }
  
  tolower(s)
}