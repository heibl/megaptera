## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-11-29)

#' @export

dbChooseRef <- function(x, spec){
  
  spec <- paste(spec, collapse = "|")
  spec.tab <- paste("spec", x@locus@sql, sep = "")
  ref <- dbReadDNA(x, tab.name = spec.tab, taxon = spec, regex = FALSE)
  ref <- del.gaps(ref)
  ref <- mafft(ref, exec = x@align.exe)
  ref
}