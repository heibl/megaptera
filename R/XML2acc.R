## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-03-24)

#' @export
#' @import XML

XML2acc <- function(x, xml, taxon){
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  logfile <- paste0("log/", gene, "-stepB.log")
  
  ## get GI with xPath
  ## -----------------
  #     gi <- xpathSApply(xml, "//GBSeqid[matches(., 'gi')]", xmlValue)
  ## the preceding line should be preferred, but xPath2 doesn't seem
  ## to be implemented in libxml2 ???
  gi <- xpathApply(xml, "//GBSeq_other-seqids", xmlToList)
  gi <- lapply(gi, unlist)
  gi <- sapply(gi, function(x) x[grep("^gi[|]", x)])
  gi <- gsub("^gi[|]", "", gi)
  
  ## get organism = taxon name
  ## -------------------------
  organism <- xpathSApply(xml, "//GBSeq_organism", xmlValue)
  if (!length(grep(" " , organism)))
    organism <- paste(organism, "sp.") # sometimes only a genus name is given in the organism field
  # organism <- gsub(" ", "_", organism) # use underscores (deprecated as of 2017-01-11)
  
  ## TOPOLOGY + LENGTH
  ## sequences > 200000 bp seem not to be included in XML!
  ## -----------------------------------------------------
  topology <- xpathSApply(xml, "//GBSeq_topology", xmlValue)
  bp <- as.numeric(xpathSApply(xml, "//GBSeq_length", xmlToList))
  
  ## get DNA
  ## -------
  dna <- xpathSApply(xml, "//GBSeq_sequence", xmlToList)
  
  ## necessary because some GBSeq elements do not have
  ## a GBSeq_sequence child:
  if (length(dna) < length(gi)){
    ids <- paste("gi", gi, sep = "|")
    ids <- paste("//GBSeq/GBSeq_other-seqids[GBSeqid = '",
                 ids, "']/../GBSeq_sequence", sep = "")
    dna2 <- lapply(ids, xpathSApply, doc = xml, fun = xmlToList)
    ids <- which(sapply(dna2, length) == 1)
    dna <- dna[ids]
    gi <- gi[ids]
    organism <- organism[ids]
    topology <- topology[ids]
    bp <- bp[ids]
  }
  if (!length(dna)) return(NULL) # 
  
  ## extracted locus from annotated genome
  ## -------------------------------------
  ## Note: extraction is triggered by a sequence
  ## length exceeding 2000 bp
  id <- which(bp > 2000)
  if (length(id) & x@locus@use.genomes){
    slog(paste("\n.. extracting '", x@locus@aliases[1], 
               "' from ", length(id), 
               " annotated genomes ..", sep = ""), file = logfile)
    if (x@locus@kind %in% c("gene", "rRNA")){
      dna[id] <- sapply(gi[id], extractLocus, xml = xml, 
                        locus = x@locus)
    } else {
      dna[id] <- sapply(gi[id], extractIGS, xml = xml, 
                        locus = x@locus)
    }
  }
  
  ## if extraction failed, NA is returned
  ## ------------------------------------
  id <- sapply(dna, is.na)
  if (any(id)){
    if (length(which(id)) == length(gi)) return(NULL)
    gi <- gi[!id]
    organism <- organism[!id]
    topology <- topology[!id]
    dna <- dna[!id]
  }
  
  ## checkpoint: was there failure to retrieve a single gene
  check <- which(dna == "")
  if (length(check)){
    warning("sequence retrieval failed for", 
            paste("\n", sql.wrap(gi[check], term = "gi", BOOL = NULL)))
  }
  
  ## transform list in data frame
  ## ----------------------------
  if (is.Linnean(taxon[1])){
    taxa <- taxon[1]
    status <- "raw"
  } else {
    taxa <- organism
    indet <- indet.strings(x@taxon@hybrids, TRUE)
    indet <- grep(indet, taxa)
    status <- rep("raw", length(taxa))
    status[indet] <- "excluded (indet)"
    if (length(indet)){
      taxa[-indet] <- strip.infraspec(taxa[-indet])
    } else {
      taxa <- strip.infraspec(taxa)
    }
  }
  seqs <- data.frame(gi = gi,
                     taxon = taxa,
                     spec_ncbi = organism,
                     status = status,
                     genom = topology,
                     npos = sapply(dna, nchar),
                     identity = NA,
                     coverage = NA,
                     dna = dna,
                     stringsAsFactors = FALSE)
  
  seqs <- unique(seqs) # yes, there are duplicated UIDs returned by eSearch!
  
  ## write into pgSQL database
  ## -------------------------
  if ( nrow(seqs) > 0 ) {
    conn <- dbconnect(x@db)
    present <- dbGetQuery(conn, paste("SELECT gi FROM", acc.tab))
    if ( nrow(present) > 0 ){
      id <- seqs$gi %in% present$gi
      slog("\n..", length(which(id)), "duplicates removed ..", file = logfile)
      seqs <- seqs[!id, ]
    }
    dbWriteTable(conn, acc.tab, seqs, row.names = FALSE, append = TRUE)
    slog("\n..", nrow(seqs), "sequences written to", acc.tab, "", file = logfile)  
    dbDisconnect(conn)
  } # end of IF-clause (l.123)
}