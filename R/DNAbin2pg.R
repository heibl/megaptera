## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2018-01-30)

#' @export

DNAbin2pg <- function(DNAbin, reliability, mode = "character"){
  
  mode <- match.arg(mode, c("character", "raw"))
  
  ## Prepare sequences
  ## -----------------
  if (is.matrix(DNAbin)) DNAbin <- as.list(DNAbin)
  if (mode == "character") DNAbin <- as.character(DNAbin)
  DNAbin <- lapply(DNAbin, paste, collapse = "")
  
  ## Prepare reliability scores
  ## --------------------------
  if (!missing(reliability)) {
    reliability <- round(reliability, 2)
    ## Reliability can be ...
    if (!is.null(dim(reliability))){
      ## ... a matrix of cell scores ...
      reliability <- apply(reliability, 1, paste, collapse = " ")
    } else {
      ## ... or a vector of columns scores.
      reliability <- paste(reliability, collapse = " ")
    }
    
  } else {
    reliability <- -1
  }
  
  
  out <- data.frame(taxon = names(DNAbin), 
                    sequence = unlist(DNAbin),
                    reliability = reliability,
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}