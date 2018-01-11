## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-11-23)

#' @export

downloadSequences <- function(x, taxon){
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  logfile <- paste0("log/", gene, "-stepB.log")
  retmax <- x@params@retmax
  dbl <- x@params@debug.level
  
  slog("\n\n.. search taxon:", taxon, "..", file = logfile, megProj = x)
  
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
  xml <- gsub(" ", "+", xml)
  
  ## get and parse results via eFetch
  ## --------------------------------
  xml <- robustXMLparse(xml, logfile = logfile)
  webEnv <- xpathSApply(xml, fun = xmlToList,
                        path = "//eSearchResult/WebEnv")
  queryKey <- xpathSApply(xml, fun = xmlToList,
                          path = "//eSearchResult/QueryKey")
  
  n <- as.numeric(xpathSApply(xml, fun = xmlValue,
                              path = "//eSearchResult/Count"))
  if (n == 0) {
    slog("\n.. no sequences available ..", file = logfile, megProj = x)
    return(NULL)
  }
  slog("\n..", n, "UIDs posted on Entrez History Server ..", 
       file = logfile, megProj = x)
  
  ## do not download more than 100 sequences if taxon
  ## is of rank species or below (e.g. Manduca sexta 
  ## as outgroup with 30000 16S sequences)
  ## -------------------------------------
  if (is.Linnean(taxon) & n > 100){
    slog("\n.. reached limit of 100 sequences per species ..", 
         file = logfile, megProj = x)
    n <- 100
  }
  
  ## start downloading only if at least one sequence
  ## needs to be downloaded
  ## -----------------------------------------------
  ## NOTE as of 2016-07-26 this XML still contains only GIs and no GBs
  uid <- unlist(xpathSApply(xml, "//eSearchResult/IdList", xmlToList))
  conn <- dbconnect(x)
  present.gi <- dbGetQuery(conn, paste("SELECT gi FROM", acc.tab))$gi
  dbDisconnect(conn)
  
  ## CASE 1: THERE ARE NO NEW SEQUENCES
  ## ----------------------------------
  missing.gi <- setdiff(uid, present.gi)
  if (!length(missing.gi)) {
    slog("\n.. all", n, "sequences already in database ..", file = logfile, megProj = x)
    return(NULL)
  }
  
  ## CASE 2: SOME SEQUENCES PRESENT, SOME NOT
  ## ----------------------------------------
  already.present <- intersect(uid, present.gi)
  if (length(already.present)){
    
    ## DIRTY HACK NOTE: The indexing below is arbitrary. The
    ## https interface cannot handle more IDs. To avoid more complicated
    ## code I simply decided to truncate after 700 IDs. If the code is run 
    ## more than once succesively more sequences will be picked up. [2016-12-05]
    if (length(missing.gi) > 700){
      missing.gi <- missing.gi[1:700]
    }
    
    slog("\n..", length(already.present), 
         "sequences already in database ..\n..",
         length(missing.gi), "sequences will be downloaded ..", 
         file = logfile, megProj = x)
    xml <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
                  "efetch.fcgi?tool=megaptera&email=heibsta@gmx.net",
                  "&db=nucleotide", 
                  "&rettype=gb&retmode=xml",
                  "&id=", paste(missing.gi, collapse = ","))
    xml <- robustXMLparse(xml, logfile = logfile)
    XML2acc(x, xml, taxon)
    
  } else {
    
    ## CASE 3: ALL SEQUENCES ARE NEW
    ## -----------------------------
    
    ## loop over sliding window
    ## ------------------------
    sw <- seq(from = 0, to = n, by = retmax)
    sw <- data.frame(from = sw, to = c(sw[-1] - 1, n))
    b <- ifelse(nrow(sw) == 1, "batch", "batches")
    slog("\n.. retrieving records in", nrow(sw), b, "..", 
         file = logfile, megProj = x)
    
    for (i in 1:nrow(sw)) {
      
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
      if (x@params@debug.level > 3){
        err_file <- paste0("log/", Sys.Date(), "_DEBUG_downloadingSequences.rda")
        save(x, xml, file = err_file)
      }
      if (is.null(xml)) next
      # saveXML(xml, "megapteraAP-DEBUG.xml"); system("open -t megapteraAP-DEBUG.xml")
      XML2acc(x, xml, taxon)
    } # end of FOR-loop
  } # end of CASE 3
}