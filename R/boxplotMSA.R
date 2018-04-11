## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-02-24)

#' @export
#' @importFrom graphics boxplot
#' @importFrom stats median

boxplotMSA <- function(x){
  
  tabs <- dbTableNames(x, x@taxon@tip.rank)
  
  dbDistDNA <- function(msa, model = "JC69"){
    
    s <- dbReadDNA(x, msa, taxon = ".*")
    if ( is.null(s) ) return(NULL)
    dist.dna(s, model = "raw", pairwise.deletion = TRUE)
  }
  obj <- lapply(tabs, dbDistDNA)
  names(obj) <- tabs
  obj <- obj[!sapply(obj, is.null)]
  obj <- obj[order(sapply(obj, median))]
  
  obj <- lapply(obj, as.vector)
  names(obj) <- gsub(paste(x@taxon@tip.rank, "_", sep = ""), "", names(obj))
  names(obj) <- gsub("_", ".", names(obj))
  out <- boxplot(obj, horizontal = TRUE, las = 1)
  invisible(out)
}