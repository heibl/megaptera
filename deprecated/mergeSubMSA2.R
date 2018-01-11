## This code is part of the megaptera package
## difference to mergeSubMSA:
## - does not write MSA into database
## - does not return subtree label but MSA

#' @export

mergeSubMSA2 <- function(megProj, subtree){
   
  ## PARAMETERS
  ## -----------
  gene <- megProj@locus@sql
  msa.tab <- paste(megProj@taxon@tip.rank, 
                   gsub("^_", "", gene), 
                   sep = "_")
  merge.exe <- megProj@merge.exe
  
  ## read a pair of subMSAs from database:
  subMSA <- lapply(subtree, dbReadDNA, x = megProj, 
                   tab.name = msa.tab, regex = FALSE, subtree = TRUE)
  
  ## merge pair of subMSAs using OPAL:
  opal.merge(subMSA, merge.exe)
  
}