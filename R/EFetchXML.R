## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-04-14)

EFetchXML <- function (gi, debug = FALSE){
  
  ## retrieve sequence information
  ## -----------------------------
  gi <- paste(gi, collapse = ",")
  url <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
             "efetch.fcgi?db=nucleotide&id=", gi, 
             "&rettype=gb&retmode=xml", sep = "")
  xml <- xmlTreeParse(url, getDTD = FALSE, useInternalNodes = TRUE)
  
  if ( debug ){
    saveXML(xml, "test2.xml"); system("open -t test2.xml")
  }
    
  ## return NA if download failed
  ## ----------------------------
  err <- xpathSApply(xml, "/GBSet/Error", xmlValue)
  if ( length(err) > 0 ){
    warning(err)
    return(NA)
  } 
  xml
}