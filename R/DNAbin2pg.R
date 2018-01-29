## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2018-01-26)

#' @export

DNAbin2pg <- function(DNAbin, reliability, mode = "character"){
  
  mode <- match.arg(mode, c("character", "raw"))
  if (mode == "character") DNAbin <- as.character(DNAbin)
  
  if (is.list(DNAbin)){
    DNAbin <- lapply(DNAbin, function(z) {dim(z) <- NULL; z})
    out <- list()
    for (i in seq_along(DNAbin)){
      if (mode == "raw") class(DNAbin[[i]]) <- NULL
      out[[i]] <- data.frame(taxon = names(DNAbin)[i], 
                             nuc = DNAbin[[i]],
                             pos = 1:length(DNAbin[[i]]), 
                             reliability = -1,
                             stringsAsFactors = FALSE)
      
    }
  }
  if (is.matrix(DNAbin)){
    out <- list(length = nrow(DNAbin))
    pos <- 1:ncol(DNAbin)
    if (missing(reliability)){
      reliability <- matrix(1, nrow = nrow(DNAbin), ncol = ncol(DNAbin))
    }
    for (i in seq_along(rownames(DNAbin))){
      out[[i]] <- data.frame(taxon = rownames(DNAbin)[i], 
                             nuc = DNAbin[i, ],
                             pos = pos, 
                             reliability = reliability[i, ],
                             stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}