## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-03-28)

#' @export

speciesConsensus <- function(megProj, spec){
  
  md5 <- spec[3]
  n <- spec[2]
  spec <- spec[1]
  
  ## PARAMETERS
  ## -----------
  gene <- megProj@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  align.exe <- megProj@align.exe
  max.bp <- megProj@params@max.bp
  min.identity <- megProj@locus@min.identity
  min.coverage <- megProj@locus@min.coverage
  logfile <- paste0("log/", gene, "stepF.log")
  
  ## Read species alignments
  ## -----------------------
  obj <- dbReadDNA(spec, x = megProj, tab.name = acc.tab, 
                   max.bp = max.bp, 
                   min.identity = min.identity, 
                   min.coverage = min.coverage)
  if (!is.matrix(obj)) obj <- mafft(obj, exec = megProj@align.exe)

  ## Species consensus sequences
  ## ---------------------------
  obj <- specCons(obj, log = logfile)
  obj <- list(obj)
  names(obj) <- spec$taxon
  class(obj) <- "DNAbin"
  
  ## Write MSA to database
  ## ---------------------
  dbWriteMSA(megProj, obj, n = n, md5 = md5)
}