## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-01-12)

#' @export

coverage <- function(DNAbin, what = "species"){
  
  if ( !inherits(DNAbin, "DNAbin") )
    stop("'DNAbin' must be of class 'DNAbin'")
  if ( !is.matrix(DNAbin) )
    stop("'DNAbin' must be alignend")
  
  what <- match.arg(what, c("species", "nucleotides"))
  
  b <- as.raw(c(136, 40, 72, 24))
  percentInformation <- function(x, b){
    length(which(x %in% b))
  }
  if ( what == "species" ){
    m <- apply(DNAbin, 1, percentInformation, b)
    m <- m/ncol(DNAbin)
  } else {
    m <- apply(DNAbin, 2, percentInformation, b)
    m <- m/nrow(DNAbin)
  }
  m
}