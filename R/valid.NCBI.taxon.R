## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-09-15)

## kingdom is not accounted for
## both accepted and synonyms return TRUE
## should be able to handle more than 1 taxon

valid.NCBI.taxon <- function(taxon){
  
  term <- paste(taxon, collapse = "+OR+")
  xml <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/", 
               "eutils/esearch.fcgi?",
               "tool=megaptera",
               "&email=heibsta@gmx.net",
               "&usehistory=n", 
               "&db=taxonomy",
               "&term=", term)
  xml <- xmlParse(getURL(xml))
  
  error <- xpathSApply(xml, fun = xmlValue,
                        path = "//PhraseNotFound")
  
  if ( length(error) == 0 ){
    out <- rep(TRUE, length(taxon))
  } else {
    out <- !taxon %in% error
  }
  names(out) <- taxon
  out
}


