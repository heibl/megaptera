## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-02-26)

#' @export

extractLocus <- function(xml, seqid, locus, kind = "gene"){
  
  debug.level <- 0 ## has be implemented in megapteraPars
  
  ## prepare 'seqid': GI/GB
  ## ----------------------
  seqid <- as.character(seqid)
  if (length(grep("^[[:digit:]]", seqid))){
    ## seqid: GI
    UID <- paste("gi", seqid, sep = "|")
    UID <- paste("GBSeq_other-seqids/GBSeqid='", UID, "'", sep = "")
  } else {
    ## seqid: GB
    UID <- paste("GBSeq_locus='", seqid, "'", sep = "")
  }
  
  ## prepare locus information
  ## -------------------------
  if (inherits(locus, "locus")){
    loc <- locus@aliases
    GBFeature_key <- locus@kind
  } else {
    loc <- locus
    GBFeature_key <- kind
  }
  
  GBlocus <- wrapSQL(loc, "GBQualifier_value", "=", "or")
  if (GBFeature_key == "gene") GBFeature_key <- c(GBFeature_key, "CDS")
  GBFeature_key <- wrapSQL(GBFeature_key, "GBFeature_key", "=", "or")
  xpath <- paste0("//GBSeq[", UID, "]",
                 "//GBFeature[", GBFeature_key, "]", # or GBFeature_key='tRNA'
                 "//GBQualifier[", GBlocus, 
                 "]/../../GBFeature_intervals/GBInterval/",
                 "GBInterval_", c("from", "to"))
  
  fromto <- data.frame(
    from = as.numeric(xpathSApply(xml, xpath[1], xmlValue)),
    to = as.numeric(xpathSApply(xml, xpath[2], xmlValue)))
  fromto <- unique(fromto) ## nesessary e.g. with gi=32480822 and rbcL because 
  ## both feature_keys 'gene' and 'CDS' are present
  
  ## organism
  ## --------
  organism <- paste("//GBSeq[", UID, "]",
                    "//GBSeq_organism", sep = "")
  organism <- xpathSApply(xml, organism, xmlValue)
  
  ## if extraction fails ...
  ## -----------------------
  if (!nrow(fromto)){
    genes <- paste("//GBSeq[", UID, "]",
                   "//GBQualifier[GBQualifier_name='gene']/GBQualifier_value", 
                   sep = "")
    genes <- sort(unique(xpathSApply(xml, genes, xmlValue)))
    notes <- paste("//GBSeq[", UID, "]",
                   "//GBQualifier[GBQualifier_name='note']/GBQualifier_value", 
                   sep = "")
    notes <- sort(unique(xpathSApply(xml, notes, xmlValue)))
    if (debug.level){
      save(seqid, locus, notes, 
           file = paste("debugMe", seqid, "extractLocus.rda", 
                        sep = "-"))
      #     saveXML(xml, file = "debugMe-extractLocus.xml")
    }
    warning("failed to extract locus='", loc[1], 
            "' from seqid=", seqid, " (", organism, ")", 
            "\n\navailable genes: ", paste(genes, collapse = ", "), sep = "")
    return(NA)
  }
  
  ## extract locus; convert to complement, if necessary
  ## --------------------------------------------------
  dna <- paste0("//GBSeq[", UID, "]",
               "//GBSeq_sequence")
  dna <- xpathSApply(xml, dna, xmlValue)
  
  cutSequence <- function(seq, pos){
    seq <- substr(seq, min(pos), max(pos))
    if ( pos[1] > pos[2] ) {
      seq <- unlist(strsplit(seq, ""))
      seq <- rev(seqinr::comp(seq))
      seq <- paste(seq, collapse = "")
    }
    seq
  }
  dna <- apply(fromto, 1, cutSequence, seq = dna)
  dna <- paste(dna, collapse = "")
  names(dna) <- gsub(" ", "_", organism) ## add name of organism
  
  return(dna)
}