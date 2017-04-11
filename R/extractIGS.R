## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-04-14)

extractIGS <- function(xml, gi, locus){
  
  GBgi <- paste("gi", gi, sep = "|")
  
  if ( inherits(locus, "locus") ){
    adj1 <- locus@adj.gene1
    adj2 <- locus@adj.gene2
  } else {
    adj1 <- locus[[1]]
    adj2 <- locus[[2]]
  }
  
  ## upstream
  ## --------
  GBlocus <- sql.wrap(adj1, regex = FALSE, BOOL = "or", 
                      term = "GBQualifier_value")
  xpath <- paste("//GBSeq[GBSeq_other-seqids/GBSeqid='", GBgi, "']",
                 "//GBFeature[GBFeature_key='gene' or GBFeature_key='CDS']", # or GBFeature_key='tRNA'
                 "//GBQualifier[", GBlocus, 
                 "]/../../GBFeature_intervals/GBInterval/",
                 "GBInterval_", c("from", "to"), 
                 sep = "")
  fromto1 <- data.frame(
    from = as.numeric(xpathSApply(xml, xpath[1], xmlValue)),
    to = as.numeric(xpathSApply(xml, xpath[2], xmlValue)))
  
#   if ( nrow(fromto1) == 0 ) stop("query failed")
  
  ## downstream
  ## ----------
  GBlocus <- sql.wrap(adj2, regex = FALSE, BOOL = "or", 
                      term = "GBQualifier_value")
  xpath <- paste("//GBSeq[GBSeq_other-seqids/GBSeqid='", GBgi, "']",
                 "//GBFeature[GBFeature_key='gene']", # or GBFeature_key='tRNA'
                 "//GBQualifier[", GBlocus, 
                 "]/../../GBFeature_intervals/GBInterval/",
                 "GBInterval_", c("from", "to"), 
                 sep = "")
  fromto2 <- data.frame(
    from = as.numeric(xpathSApply(xml, xpath[1], xmlValue)),
    to = as.numeric(xpathSApply(xml, xpath[2], xmlValue)))
  
  ## handle failure to match adj.gene1 or adj.gene2
  ## ----------------------------------------------
  if ( nrow(fromto1) == 0 | nrow(fromto2) == 0 ) {
    genes <- paste("//GBSeq[GBSeq_other-seqids/GBSeqid='", GBgi, "']",
                   "//GBQualifier[GBQualifier_name='gene']/GBQualifier_value", 
                   sep = "")
    genes <- sort(unique(xpathSApply(xml, genes, xmlValue)))
    notes <- paste("//GBSeq[GBSeq_other-seqids/GBSeqid='", GBgi, "']",
                   "//GBQualifier[GBQualifier_name='note']/GBQualifier_value", 
                   sep = "")
    notes <- sort(unique(xpathSApply(xml, notes, xmlValue)))
    save(gi, adj1, adj2, notes, 
         file = paste("debugMe", gi, "extractLocus.rda", 
                      sep = "-"))
    #     saveXML(xml, file = "debugMe-extractLocus.xml")
    warning("failed to extract locus='", adj1[1], "-", adj2[1], 
            " intergenic spacer' from gi=", gi, 
            "\n\navailable genes: ", 
            paste(genes, collapse = ", "), sep = "")
    return(NA)
  }
  #   cat("\n")
  #   print(fromto1)
  #   cat("\n")
  #   print(fromto2)
  
  ## tRNAs can be present more than one time
  ## ---------------------------------------
  v <- vector()
  for ( i in rowMeans(fromto1) ){
    for ( j in rowMeans(fromto2) )
      v <- c(v, abs(i - j))
  }
  v <- matrix(v, nrow = nrow(fromto1), byrow = TRUE)
  id <- arrayInd(which.min(v), 
                 .dim = c(nrow(fromto1), nrow(fromto2)))
  fromto1 <- fromto1[id[, 1], ]
  fromto2 <- fromto2[id[, 2], ]
  
  ## detect order of adjacent genes
  if ( max(fromto1) < min(fromto2) ){
    fromto <- cbind(max(fromto1), min(fromto2))
  } else {
    fromto <- cbind(max(fromto2), min(fromto1))
  }
  
  ## detect if spacer is on opposing strand
  if ( fromto1[1, 1] > fromto1[1, 2] ) {
    fromto <- fromto[, 2:1, drop = FALSE]
  }
  
  print(fromto)
  
  ## detect some caveats
  ## -------------------
  if (  nrow(fromto) == 0 ){
    loci <- paste("//GBSeq[GBSeq_other-seqids/GBSeqid='", GBgi, "']",
                  "//GBQualifier[GBQualifier_name='gene']/GBQualifier_value", 
                  sep = "")
    loci <- sort(unique(xpathSApply(xml, loci, xmlValue)))
    stop("failed to extract locus='", locus[1], "' from gi=", gi, 
         "\navailable loci: ", paste(loci, collapse = ", "), sep = "")
  }
  
  ## extract locus; convert to complement, if necessary
  ## --------------------------------------------------
  dna <- paste("//GBSeq[GBSeq_other-seqids/GBSeqid='", GBgi, "']",
               "//GBSeq_sequence", 
               sep = "")
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
  
  ## add name of organism
  organism <- paste("//GBSeq[GBSeq_other-seqids/GBSeqid='", GBgi, "']",
                    "//GBSeq_organism", sep = "")
  organism <- xpathSApply(xml, organism, xmlValue)
  names(dna) <- gsub(" ", "_", organism)
  
  return(dna)
}