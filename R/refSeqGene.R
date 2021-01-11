## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-11-24)

#' @import DBI
#' @importFrom RCurl url.exists
#' @export

refSeqGene <- function(locus){
  
  ## CHECKS
  ## ------
  if (!url.exists("https://eutils.ncbi.nlm.nih.gov"))
    stop("internet connection required")
  
  ## Get XML
  ## ----------------------------------
  xml <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/", 
                "eutils/esearch.fcgi?",
                "tool=megaptera",
                "&email=heibsta@gmx.net",
                "&db=nucleotide",
                "&term=(RefSeqGene[keyword]+", 
                "AND+", locus, "[Gene]")
  xml <- gsub(" ", "+", xml)
  xml <- robustXMLparse(xml)
  xml <- unlist(xpathSApply(xml, "//eSearchResult/IdList", xmlToList))
  
  if (is.null(xml)) stop("no reference sequences found")
  if (length(xml) > 1){
    def <- paste(xml, collapse = ",")
    def <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
                  "efetch.fcgi?db=nucleotide", 
                  "&id=", def, 
                  "&rettype=gb", 
                  "&retmode=xml")
    def <- xmlTreeParse(getURL(def), 
                        getDTD = FALSE, 
                        useInternalNodes = TRUE)
    def <- xpathSApply(def, "//GBSeq_definition", xmlValue)
    message("Choices:")
    for (i in 1:length(def)){
      message(paste0("(", i, ") ", def[i]))
    }
    id <- as.numeric(readline())
    xml <- xml[id]
  }
  
  xml <- paste(xml, collapse = ",")
  xml <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
                "efetch.fcgi?db=nucleotide", 
                "&id=", xml, 
                "&rettype=gb", 
                "&retmode=xml")
  xml <- xmlTreeParse(getURL(xml), 
                      getDTD = FALSE, 
                      useInternalNodes = TRUE)
  message("Using: ", xpathSApply(xml, "//GBSeq_definition", xmlValue))
  
  # saveXML(xml, "test2.xml"); system("open -t test2.xml")
  
  ## Parse XML (see APPENDIX A)
  ## --------------------------
  organism <- xpathSApply(xml, "//GBSeq_organism", xmlValue)
  acc <- xpathSApply(xml, "//GBSeq_primary-accession", xmlValue)
  
  xpath <- paste0("//GBSeq",
                  "//GBFeature[GBFeature_key='gene']", # or 
                  "/GBFeature_quals/GBQualifier[GBQualifier_value='", locus, "']",
                  "/../../GBFeature_intervals/GBInterval/",
                  "GBInterval_", c("from", "to"))
  fromto <- data.frame(
    from = as.numeric(xpathSApply(xml, xpath[1], xmlValue)),
    to = as.numeric(xpathSApply(xml, xpath[2], xmlValue)))
  
  dna <- xpathSApply(xml, "//GBSeq//GBSeq_sequence", xmlValue)
  dna <- substr(dna, fromto$from, fromto$to)
  
  ## Create DNAbin
  obj <- strsplit(dna, "")
  names(obj) <- gsub(" ", "_", paste(organism, acc))
  obj <- as.DNAbin(obj)
  
  obj
}
# APPENDIX A
# <GBFeature>
#   <GBFeature_key>gene</GBFeature_key>
#   <GBFeature_location>5001..47645</GBFeature_location>
#   <GBFeature_intervals>
#     <GBInterval>
#       <GBInterval_from>5001</GBInterval_from>
#       <GBInterval_to>47645</GBInterval_to>
#       <GBInterval_accession>NG_011793.1</GBInterval_accession>
#     </GBInterval>
#   </GBFeature_intervals>
#   <GBFeature_quals>
#   <GBQualifier>
#   <GBQualifier_name>gene</GBQualifier_name>
#   <GBQualifier_value>APOB</GBQualifier_value>
#   </GBQualifier>
#   <GBQualifier>
#   <GBQualifier_name>gene_synonym</GBQualifier_name>
#   <GBQualifier_value>apoB-100; apoB-48; FCHL2; FLDB; LDLCQ4</GBQualifier_value>
#   </GBQualifier>
#   <GBQualifier>
#   <GBQualifier_name>note</GBQualifier_name>
#   <GBQualifier_value>apolipoprotein B</GBQualifier_value>
#   </GBQualifier>
#   <GBQualifier>
#   <GBQualifier_name>db_xref</GBQualifier_name>
#   <GBQualifier_value>GeneID:338</GBQualifier_value>
#   </GBQualifier>
#   <GBQualifier>
#   <GBQualifier_name>db_xref</GBQualifier_name>
#   <GBQualifier_value>HGNC:HGNC:603</GBQualifier_value>
#   </GBQualifier>
#   <GBQualifier>
#   <GBQualifier_name>db_xref</GBQualifier_name>
#   <GBQualifier_value>MIM:107730</GBQualifier_value>
#   </GBQualifier>
#   </GBFeature_quals>
#   </GBFeature>