## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-08-16)

mergeSubMSA <- function(x, subtree){
   
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  msa.tab <- paste(x@taxon@tip.rank, 
                   gsub("^_", "", gene), 
                   sep = "_")
  merge.exe <- x@merge.exe
  subMSA <- lapply(subtree, dbReadDNA, x = x, 
                   tab.name = msa.tab, subtree = TRUE)
  obj <- opal.merge(subMSA, merge.exe)
  
  subtree <- sort(subtree)
  subtree[2] <- gsub("st-", "", subtree[2])
  subtree <- paste(subtree, collapse = "-")
  dbWriteMSA(x, obj, status = "merged", subtree = subtree)
  subtree
}