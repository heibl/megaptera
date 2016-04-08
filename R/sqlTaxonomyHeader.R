## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-07-24)

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
  names(tax) <- gsub("family", "fam", names(tax))
  names(tax) <- gsub("order", "ord", names(tax))
  names(tax) <- gsub("class", "class", names(tax))
  names(tax) <- gsub("[.]", "_", names(tax))
  
  if ( mode == "reference.rank" ){
    tax <- names(tax)
  }
  tax
}
