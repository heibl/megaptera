## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-09-15)

downloadSequences <- function(x, taxon){
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  logfile <- paste(gene, "stepB.log", sep = "-")
  retmax <- x@params@retmax
  
  slog("\n\n.. search taxon:", taxon, "..", file = logfile)
  
  ## post UIDs on Entrez History Server
  ## ----------------------------------
  xml <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/", 
               "eutils/esearch.fcgi?",
               "tool=megaptera",
               "&email=heibsta@gmx.net",
               "&usehistory=y", 
               "&retmax=99999",
               "&db=nucleotide")
  xml <- paste0(xml, "&term=", 
               term(taxon, x@taxon@kingdom, x@locus))
  
  ## get and parse results via eFetch
  ## --------------------------------
  xml <- robustXMLparse(xml, logfile = logfile)
  webEnv <- xpathSApply(xml, fun = xmlToList,
                        path = "//eSearchResult/WebEnv")
  queryKey <- xpathSApply(xml, fun = xmlToList,
                          path = "//eSearchResult/QueryKey")
  
  n <- as.numeric(xpathSApply(xml, fun = xmlValue,
                   path = "//eSearchResult/Count"))
  if ( n == 0 ) {
    slog("\n.. no sequences available ..", file = logfile)
    return(NULL)
  }
  
  ## start downloading only if at least one sequence
  ## needs to be downloaded
  ## -----------------------------------------------
  ## NOTE as of 2016-07-26 this XML still contains only GIs and no GBs
  uid <- unlist(xpathSApply(xml, "//eSearchResult/IdList", xmlToList))
  conn <- dbconnect(x)
  present.gi <- dbGetQuery(conn, paste("SELECT gi FROM", acc.tab))$gi
  dbDisconnect(conn)
  if ( length(setdiff(uid, present.gi)) == 0 ) {
    slog("\n.. all", n, "sequences already in database ..", file = logfile)
    return(NULL)
  }
  
  ## loop over sliding window
  ## ------------------------
  sw <- seq(from = 0, to = n, by = retmax)
  sw <- data.frame(from = sw, to = c(sw[-1] - 1, n))
  slog("\n.. posting", n, "UIDs on Entrez History Server ..", 
       file = logfile)
  b <- ifelse(nrow(sw) == 1, "batch", "batches")
  slog("\n.. retrieving full records in", nrow(sw), b, "..", 
       file = logfile)
  
  for ( i in 1:nrow(sw) ) {
    
    ## get XML with full records
    ## -------------------------
    xml <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
                 "efetch.fcgi?tool=megaptera&email=heibsta@gmx.net",
                 "&db=nucleotide&query_key=", queryKey, 
                 "&WebEnv=", webEnv,
                 "&rettype=gb&retmode=xml",
                 "&retstart=", sw$from[i], 
                 "&retmax=", retmax)
    
    ## parse XML: this step is error-prone and therefore
    ## embedded into try()
    ## -------------------
    xml <- robustXMLparse(xml, logfile = logfile)
    if ( is.null(xml) ) next
    
    # saveXML(xml, "megapteraAP-DEBUG.xml"); system("open -t megapteraAP-DEBUG.xml")
     
    
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
    if ( length(grep(" " , organism)) == 0 )
      organism <- paste(organism, "sp.") # sometimes only a genus name is given in the organism field
    organism <- gsub(" ", "_", organism) # use underscores
    
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
    if ( length(dna) < length(gi) ){
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
    if ( length(dna) == 0 ) next # 
    
    ## extracted locus from annotated genome
    ## -------------------------------------
    ## Note: extraction is triggered by a sequence
    ## length exceeding 2000 bp
    id <- which( bp > 2000 )
    if ( length(id) > 0 & x@locus@use.genomes ){
      slog(paste("\n.. extracting '", x@locus@aliases[1], 
                 "' from ", length(id), 
                 " annotated genomes ..", sep = ""), file = logfile)
      if ( x@locus@kind %in% c("gene", "rRNA") ){
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
    if ( any(id) ){
      if ( length(which(id)) == length(gi) ) return(NULL)
      gi <- gi[!id]
      organism <- organism[!id]
      topology <- topology[!id]
      dna <- dna[!id]
    }
    
    ## checkpoint: was there failure to retrieve a single gene
    check <- which(dna == "")
    if ( length(check) > 0 ){
      warning("sequence retrieval in batch", i, "failed for", 
              paste("\n", sql.wrap(gi[check], term = "gi", BOOL = NULL)))
    }
    
    ## transform list in data frame
    ## ----------------------------
    if ( is.Linnean(taxon[1]) ){
      taxa <- gsub(" ", "_", taxon[1])
    } else {
      taxa <- organism
    }
#     cat("\ngi      ", length(gi))
#     cat("\ntaxon   ", length(taxon))
#     cat("\norganism", length(organism))
#     cat("\ngenom   ", length(topology))
#     cat("\ndna     ", length(dna))
    seqs <- data.frame(gi = gi,
                       taxon = taxa,
                       spec_ncbi = organism,
                       status = "raw",
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
    } # end of IF-clause (l.187)
  } # end of FOR-loop
}