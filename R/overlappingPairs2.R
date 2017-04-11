

overlappingPairs2 <- function(seqs){
  
  n <- lapply(seqs, names)
  
  obj <- matrix(0, nrow = length(seqs), ncol = length(seqs))
  for ( i in seq_along(n) ){
    id <- which(sapply(n, function(y, z) any(y %in% z), y = n[[i]]))
    id <- setdiff(id, i)
    for ( j in id ){
      obj[i, j] <- obj[j, i] <- 1
    }
  }
  o <- obj
  on.board <- vector()
  repeat {
    # for ( h in 1:7){
    cs <- colSums(obj)
    if ( all(cs == 0) ){
      id <- which(!test)
      if ( length(id) > 1 ) stop("debug me")
      a <- which(o[id[1], ] == 1)[1]
      on.board <- c(on.board, id[1], a)
      break
    }
    id <- which(cs == 1)
    if ( length(id) > 0 ){
      a <- which(obj[id[1], ] == 1)[1]
      on.board <- c(on.board, id[1], a)
      obj[c(id[1], a), ] <- obj[, c(id[1], a)] <- 0
    } else {
      id <- which(cs > 1)
      a <- which(obj[id[1], ] == 1)[1]
      on.board <- c(on.board, id[1], a)
      obj[c(id[1], a), ] <- obj[, c(id[1], a)] <- 0
    }
    test <- 1:length(seqs) %in% on.board
    cat("\n", length(which(test)))
    if ( all(test) ) break
    }
  
  op <- data.frame(a = on.board[c(TRUE, FALSE)],
                   b = on.board[c(FALSE, TRUE)])
  op <- apply(op, 1, as.list)
  lapply(op, unlist)
}

overlappingPairs <- function(x){
  
  obj <- vector(mode = "list")
  for ( i in seq_along(x) ){
    id <- which(sapply(x, function(y, z) any(y %in% z), y = x[[i]]))
    id <- setdiff(id, i)
    for ( j in id ){
      obj <- c(obj, list(sort(c(i, j))))
    }
  }
  unique(obj)
}



# min(setdiff(id, unlist(i)))

# p <- list(c(i, min(setdiff(id, unlist(i)))))



