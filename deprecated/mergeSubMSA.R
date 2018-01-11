## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-02-17)

#' @export

mergeSubMSA <- function(x, subtree){
   
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  msa.tab <- paste(x@taxon@tip.rank, 
                   gsub("^_", "", gene), 
                   sep = "_")
  merge.exe <- x@merge.exe
  
  ## read a pair of subMSAs from database:
  subMSA <- lapply(subtree, dbReadDNA, x = x, tab.name = msa.tab, 
                   regex = FALSE, subtree = TRUE)
  
  ## merge pair of subMSAs using OPAL:
  obj <- opal.merge(subMSA, merge.exe)
  
  ## create new subtree label:
  subtree <- sort(subtree)
  subtree[2] <- gsub("st-", "", subtree[2])
  subtree <- paste(subtree, collapse = "-")
  
  ## write new MSA to database:
  dbWriteMSA(x, obj, status = "merged", subtree = subtree)
  
  ## return new subtree label:
  subtree
}