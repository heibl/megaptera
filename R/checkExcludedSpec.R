## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-01-11)

checkExcludedSpec <- function(x, latex){
  
  ## get locus table with only 'blocks' columns:
  lcs <- dbReadLocus(x)
  colnames(lcs) <- gsub("sel_", "", colnames(lcs))

  ## species that contain excluded sequences/loci:
  id <- which(lcs == "excluded (by user)", arr.ind = TRUE)
  lcs <- lcs[sort(id[, "row"]), ]
  
  ## format table:
  colnames(lcs) <- gsub("_blocks", "", colnames(lcs))
  lcs[is.na(lcs)] <- 0
  lcs[lcs == "selected (block-1)"] <- 1
  lcs[lcs == "excluded (by user)"] <- "x"
  z <- lcs
  z[z == "x"] <- 0
  z <- rowSums(apply(z, c(1, 2), as.numeric))
  id <- which(lcs == "x", arr.ind = TRUE)
  lcs <- lcs[, unique(id[, "col"])]
  lcs <- data.frame(lcs, Remaining_loci = z)
  
  if ( !missing(latex) ){
    n <- ncol(lcs) - 1
    header <- c("Species", gsub("^X_", "", colnames(lcs)))
    header <- gsub("_", " ", header)
    header <- paste(header, collapse  = "&")
    header <- paste(header, "\\\\", sep = "")
    spec <- gsub("_", " ", rownames(lcs))
    spec <- paste("\\textit{", spec, "}", sep = "")
    lcs <- apply(lcs, 1, paste, collapse = "&")
    lcs <- paste(spec, "&", lcs, "\\\\", sep = "")
    tex <- paste("\\begin{tabular}{l", paste(rep("c", n), collapse = " "), "r}")
    tex <- c(tex, "\\hline", header, "\\hline", lcs, "\\hline", "\\end{tabular}")
    write(tex, latex)
  } else {
    return(lcs)
  }
}