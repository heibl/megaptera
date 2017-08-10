## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-08-08)

#' @export

EFetchXML <- function (gi, debug = FALSE){
  
  ## retrieve sequence information
  ## -----------------------------
  gi <- paste(gi, collapse = ",")
  xml <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
             "efetch.fcgi?db=nucleotide", 
             "&id=", gi, 
             "&rettype=gb", 
             "&retmode=xml")
  xml <- xmlTreeParse(getURL(xml), 
                      getDTD = FALSE, 
                      useInternalNodes = TRUE)
  
  if (debug){
    saveXML(xml, "test2.xml"); system("open -t test2.xml")
  }
    
  ## return NA if download failed
  ## ----------------------------
  err <- xpathSApply(xml, "/GBSet/Error", xmlValue)
  if (length(err)){
    warning(err)
    return(NA)
  } 
  xml
}