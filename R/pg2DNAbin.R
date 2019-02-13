## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2019-02-13)

#' @export

pg2DNAbin <- function(pg, reliability = 0){
  
  ## Parse DNA sequences
  ## -------------------
  obj <- strsplit(pg$sequence, split = "")
  names(obj) <- pg$taxon
  obj <- as.DNAbin(obj)
  if (length(unique(sapply(obj, length))) == 1){
    obj <- as.matrix(obj)
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
      rel_score <- rel_score[, rel_score >= reliability]
    }
  }
  
  
  
  attr(obj, "cs") <- rel_score
  obj
}
