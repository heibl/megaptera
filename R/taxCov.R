taxCov <- function(obj, ig, og){
  taxon_set <- unique(unlist(lapply(obj, rownames)))
  ii <- length(intersect(ig, taxon_set))
  oo <- length(intersect(og, taxon_set))
  rb <- red $ bold; mb <- magenta $ bold
  cat(silver("Taxonomic coverage:\n"))
  if (ii == 0) { col_switch <- rb } else { col_switch <- mb }
  cat(silver(" > " %+% col_switch(ii) 
             %+% " ingroup taxa (" 
             %+% magenta(round(ii / length(ig) * 100, 2)) 
             %+% "%)\n"))
  
  if (oo == 0) { col_switch <- rb } else { col_switch <- mb }
  cat(silver(" > " %+% col_switch(oo) 
             %+% " ingroup taxa (" 
             %+% magenta(round(oo / length(og) * 100, 2)) 
             %+% "%)\n"))
  cat(silver(paste(rep("-", 60), collapse = "") %+% "\n"))
  invisible(taxon_set)
}
