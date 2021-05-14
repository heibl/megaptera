

outgroup <- function(megProj, subset, sep = "_"){
  
  sep_set <- c(" ", "_")
  sep <- match.arg(sep, sep_set) ## accepted seperators
  og <- dbReadTaxonomy(megProj, subset = subset)
  og_set <- megProj@taxon@outgroup
  if (!missing(subset)){
    og_set <- intersect(og_set, og$taxon)
  }
  og <- lapply(og_set, taxdumpChildren,
               tax = og, tip.rank = megProj@taxon@tip.rank)
  og <- do.call(rbind, og)$taxon
  gsub(sep_set[!sep_set %in% sep],
       sep_set[sep_set %in% sep],
       og)
}
