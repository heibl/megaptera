## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2018-02-26)

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
  
  ## Create a subset of columns with a reliability score
  ## equal or greater 'reliability'
  ## ------------------------------
  if (reliability > 0){
    if (pg$reliability[1] != "-1"){
      rel_score <- strsplit(pg$reliability, split = " ")
      rel_score <- lapply(rel_score, as.numeric)
      rel_score <- do.call(rbind, rel_score)
      rel_score <- colMeans(rel_score)
      obj <- obj[, rel_score >= reliability]
    }
  }
  
  obj
}
