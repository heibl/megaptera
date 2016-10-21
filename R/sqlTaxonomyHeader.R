## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-09-28)

sqlTaxonomyHeader <- function(tax){
  
  mode <- "table"
  if ( is.vector(tax) ){
    mode <- "reference.rank"
    names(tax) <- tax
  }
  ## enforce standard column names
  ## -----------------------------
  names(tax) <- tolower(names(tax))
  names(tax) <- gsub("species", "spec", names(tax))
  names(tax) <- gsub("genus|genera", "gen", names(tax))
  names(tax) <- gsub("familia|family", "fam", names(tax))
  names(tax) <- gsub("ordo|order", "ord", names(tax))
  names(tax) <- gsub("classis|class", "class", names(tax))
  names(tax) <- gsub("[.]", "_", names(tax))
  
  if ( mode == "reference.rank" ){
    tax <- names(tax)
  }
  tax
}
