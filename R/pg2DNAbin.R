## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2019-05-04)

#' @export

pg2DNAbin <- function(pg, reliability = 0){
  
  ## CHECKS
  ## ------
  if (!is.data.frame(pg)) stop("'pg' must be of class 'data.frame'")
  id <- c("taxon", "sequence", "reliability") %in% names(pg)
  if (!all(id)) stop("'pg' is malformatted")
  
  ## Parse DNA sequences
  ## -------------------
  obj <- strsplit(pg$sequence, split = "")
  names(obj) <- pg$taxon
  obj <- as.DNAbin(obj)
  if (length(unique(sapply(obj, length))) == 1){
    obj <- as.matrix(obj)
  } else {
    # ## If sequences are not aligned, they cannot
    # ## have reliability scores assigned to them
    # return(obj)
  }
  
  ## If reliability scores are not available,
  ## stop here and return sequences
  ## ------------------------------
  test <- sapply(pg$reliability, nchar) == 1
  if (any(is.na(test))){
    return(obj) ## reliability is NA
  }
  if (any(test)){
    return(obj) ## reliability is a single number
  }
  
  ## Parse reliability score
  ## -----------------------
  rel_score <- strsplit(pg$reliability, split = " ")
  rel_score <- lapply(rel_score, as.numeric)
  rel_score <- do.call(rbind, rel_score)
  rel_score <- colMeans(rel_score, na.rm = TRUE)
  
  ## Check if size of MSA and CS match
  ## ---------------------------------
  if (ncol(obj) != length(rel_score))
    stop("debug me!")
  
  ## Create a subset of columns with a reliability score
  ## equal or greater 'reliability'
  ## ------------------------------
  if (reliability > 0){
    if (pg$reliability[1] != "-1"){
      obj <- obj[, rel_score >= reliability]
      rel_score <- rel_score[rel_score >= reliability]
    }
  }
  attr(obj, "cs") <- rel_score
  
  
  obj
}
