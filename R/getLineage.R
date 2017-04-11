## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-11-03)

#' @export
#' @import XML

getLineage <- function(xml, spec, kingdom){
  
  # this xPath expression avoids the 'double kingdom error' CH 2014-10-23
  xPath <- paste("/ncbi:TaxaSet/Taxon[ScientificName='", spec, 
                 "']/LineageEx[Taxon[Rank='kingdom']/ScientificName='", 
                 kingdom,"']/Taxon/", sep = "")
  out<- data.frame(name = xpathSApply(xml, paste(xPath,  "ScientificName", sep = ""), xmlValue, 
                                      namespaces =  xmlNamespaceDefinitions(xmlRoot(xml), simplify = TRUE)), 
                   rank = xpathSApply(xml, paste(xPath,  "Rank", sep = ""), xmlValue, 
                                      namespaces =  xmlNamespaceDefinitions(xmlRoot(xml), simplify = TRUE)),
                   stringsAsFactors = FALSE)
  if ( nrow(out) == 0 ){
    # eg. Polypodium hydriforme (Hydrozoa) when searching Polypodium
    return(NULL)
  } else {
    return(rbind(out, c(spec, "species")))
  }
}
