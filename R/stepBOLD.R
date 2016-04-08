## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-12-18)

stepBOLD <- function(x, overwrite = TRUE){
  
  spec <- dbReadTaxonomy(x)$spec
  spec <- gsub("_", " ", spec)
  
  ## sliding window breaks specis names vector in
  ## batches of 100
  sw <- seq(from = 1, to = length(spec), by = 100)
  sw <- data.frame(from = sw, to = c(sw[-1] -1, length(spec)))
  sw <- apply(sw, 1, function(z, spec) spec[z[1]:z[2]], spec = spec)
  
  ## wrap bold_seq to get DNAbin
  formatBOLD <- function(spec){
    b <- bold_seq(spec)
    b <- do.call(rbind, b)
    bb <- paste(b[, "name"], b[, "id"])
    bb <- gsub(" ", "_", bb)
    b <- as.list(b[, "sequence"])
    b <- lapply(b, strsplit, split = "")
    b <- lapply(b, unlist)
    b <- lapply(b, tolower)
    names(b) <- bb
    as.DNAbin(b)
  }
  
  ## download and format sequences
  b <- lapply(sw, formatBOLD)
  b <- do.call(c, b)
  
  ## manage duplicates name + id combinations
  ## where do they come from???
  d <- duplicated(names(b))
  if ( any(d) ){
    names(b)[d] <- paste(names(b)[d], 1:length(which(d)), sep = "-")
    # a <- b[names(b) == "Locusta_migratoria_CYTC5284-12"]
    # a <- mafft(a, path = x@align.exe)
    # rownames(a) <- paste(rownames(a), 1:nrow(a), sep = "_")
    # write.nex(a, "aaa.nex")
  }
  stepBX(x, b, tag = "bold", overwrite = overwrite)
}