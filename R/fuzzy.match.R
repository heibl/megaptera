## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-02-27)

fuzzy.match <- function(x, table){
  
  ## one character is differrent
  x1 <- strsplit(x, "")
  x1 <- rep(x1, length(x1[[1]]))
  for ( i in seq_along(x1) ){
    x1[[i]][i] <- "[[:alpha:]]"
  }
  x1 <- lapply(x1, paste, collapse = "")
  x1 <- lapply(x1, function(x) paste("^", x, "$", sep = ""))
  x1 <- paste(unlist(x1), collapse = "|")
  
  ## one character is erroneously added
  x2 <- strsplit(x, "")
  x2 <- rep(x2, length(x2[[1]]))
  for ( i in seq_along(x2) ){
    x2[[i]][i] <- ""
  }
  x2 <- lapply(x2, paste, collapse = "")
  x2 <- lapply(x2, function(y) paste("^", y, "$", sep = ""))
  x2 <- paste(unlist(x2), collapse = "|")
  
  ## one character is erroneously added
  x3 <- strsplit(x, "")
  x3 <- rep(x3, length(x3[[1]]))
  for ( i in seq_along(x3) ){
    x3[[i]][i] <- paste(x3[[i]][i], "[[:alpha:]]", sep = "")
  }
  x3 <- lapply(x3, paste, collapse = "")
  x3 <- lapply(x3, function(x) paste("^", x, "$", sep = ""))
  x3 <- paste(unlist(x3), collapse = "|")
  
  ## one character is erroneously added
  x4 <- strsplit(x, "")
  x4 <- rep(x4, length(x4[[1]]) - 2)
  for ( i in 2:(length(x4) + 1) ){
    z <- x4[[i - 1]][i]
    x4[[i - 1]][i] <- x4[[i - 1]][i + 1]
    x4[[i - 1]][i + 1] <- z
  }
  x4 <- lapply(x4, paste, collapse = "")
  x4 <- lapply(x4, function(x) paste("^", x, "$", sep = ""))
  x4 <- paste(unlist(x4), collapse = "|")
  
  grep(paste(x1, x2, x3, x4, sep = "|"), table)
}

# gen <- c("Oxalis", "Protea", "Crocus", "Aster")
# 
# test <- c("Proteas", "Prolea", "Proteo", "Proea",
#           "Prostea", "Protae")
# sapply(test, fuzzy.match, table = gen)

