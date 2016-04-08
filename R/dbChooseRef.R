dbChooseRef <- function(x, spec){
  
  spec <- paste(spec, collapse = "|")
  spec.tab <- paste("spec", x@locus@sql, sep = "")
  ref <- dbReadDNA(x, spec.tab, spec)
  ref <- del.gaps(ref)
  ref <- mafft(ref, path = x@align.exe)
  ref
}