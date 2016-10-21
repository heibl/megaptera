## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-09-15)

## to do:
## must be adapted for plants e.g. lines 15, 104
## extend outgroup to genus if it is a single species

ncbiGenome <- function(x, organelle, n = 5){
  
  ## set organelle
  ## -------------
  organelle <- match.arg(organelle, c("mitochondrion", "chloroplast"))
  slog("\n.. organelle     :", organelle)
  ## tag can either be 'mitochondrion' or 'mitochondrial DNA'
  organelle <- gsub("drion", "dri*", organelle)
  
  ## set ingroup root
  ## ----------------
  if ( length(x@taxon@ingroup) == 1 ){
    ingroup.root <- unlist(x@taxon@ingroup)
  } else {
    ingroup.root <- dbReadTaxonomy(x, tag = "ingroup")
    ingroup.root <- findRoot(ingroup.root)
  }
  slog("\n.. ingroup root  :", ingroup.root)
  
  ## set outgroup root
  ## -----------------
  if ( length(x@taxon@outgroup) == 1 ){
    outgroup.root <- unlist(x@taxon@outgroup)
  } else {
    outgroup.root <- dbReadTaxonomy(x, tag = "outgroup")
    outgroup.root <- findRoot(outgroup.root)
  }
  slog("\n.. outgroup root :", outgroup.root)
  outgroup.root <- gsub(" ", "+", outgroup.root)
  
  ## query ENTREZ and save history on server
  ## ---------------------------------------
  xml <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/", 
                "eutils/esearch.fcgi?",
                "tool=megaptera",
                "&email=heibsta@gmx.net",
                "&usehistory=y",
                "&retmax=9999",
                "&db=nucleotide",
                "&term=(", ingroup.root, 
                "[Organism]+AND+", organelle,
                "[Title]+AND+complete+genome[Title])",
                "+OR+(", outgroup.root, 
                "[Organism]+AND+", organelle,
                "[Title]+AND+complete+genome[Title])"
  )
  # cat(xml)
  
  ## get and parse results via eFetch from history server
  ## ----------------------------------------------------
  xml <- robustXMLparse(xml, logfile = "")
  webEnv <- xpathSApply(xml, fun = xmlToList,
                        path = "//eSearchResult/WebEnv")
  queryKey <- xpathSApply(xml, fun = xmlToList,
                          path = "//eSearchResult/QueryKey")
  nn <- as.numeric(xpathSApply(xml, fun = xmlValue,
                               path = "//eSearchResult/Count"))
  if ( nn == 0 ){
    stop("no ", organelle, " genomes available for ", x@taxon@ingroup)
  }
  
  ## loop over sliding window
  ## ------------------------
  retmax = 50
  sw <- seq(from = 0, to = nn, by = retmax)
  sw <- data.frame(from = sw, to = c(sw[-1] - 1, nn))
  slog("\n.. posting", nn, "UIDs on Entrez History Server ..")
  b <- ifelse(nrow(sw) == 1, "batch", "batches")
  slog("\n.. retrieving full records in", nrow(sw), b, "..\n")
  
  y <- data.frame()
  for ( i in 1:nrow(sw) ) {
    
    ## get XML with full records
    ## -------------------------
    xml <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
                  "esummary.fcgi?tool=megaptera&email=heibsta@gmx.net",
                  "&db=nucleotide&query_key=", queryKey, 
                  "&WebEnv=", webEnv,
                  "&retmode=xml",
                  "&retstart=", sw$from[i], 
                  "&retmax=", retmax)
    
    ## parse XML: this step is error-prone and therefore
    ## embedded into try()
    ## -------------------
    xml <- robustXMLparse(xml)
    # saveXML(xml, "aaa.xml"); system("open -t aaa.xml")
    if ( is.null(xml) ) next
    
    gb <- xpathApply(xml, "//Item[@Name = 'Caption']", xmlToList)
    gb <- sapply(gb, function(x) x$text)
    gi <- xpathApply(xml, "//Item[@Name = 'Gi']", xmlToList)
    gi <- sapply(gi, function(x) x$text)
    ti <- xpathApply(xml, "//Item[@Name = 'Title']", xmlToList)
    ti <- sapply(ti, function(x) x$text)
    ti <- gsub(" mitochondri.+complete genome", "", ti)
    ti <- data.frame(taxon = ti, gb = gb, gi = gi,
                     stringsAsFactors = FALSE)
    y <- rbind(y, ti)
  }
  
  ## delete undetermined accessions
  ## false decision: "Eucryptorrhynchus chinensis voucher ECHIN20150110"
  ## ----------------------------
  indet <- grep(paste(indet.strings(), collapse = "|"), y$taxon)
  if ( length(indet) > 0 ){
    y <- y[-indet, ]
  }
  y$taxon <- strip.infraspec(y$taxon)
  y$taxon <- gsub(" ", "_", y$taxon)
  y <- y[!duplicated(y$taxon), ]
  
  ## create a phylogenetically balanced sample
  ## of genomes in case the number of available
  ## genomes exceeds n (the number of desired genomes)
  ## -------------------------------------------------
  if ( nrow(y) > n ){
    tr <- ncbiTaxonomy(as.list(y$taxon), 
                       kingdom = x@taxon@kingdom,
                       megapteraProj = x)
    tr <- tax2tree(tr)
    tr <- compute.brlen(tr)
    # tr <- cophenetic.phylo(tr)
    branching.order <- names(sort(branching.times(tr), decreasing = TRUE))
    branching.order <- as.numeric(branching.order)
    dd <- lapply(branching.order, descendants, phy = tr, type = "d")
    nl <- sapply(dd, length)
    nl[-1] <- nl[-1] - 1
    nl <- data.frame(node = branching.order,
                     number.lineages = cumsum(nl))
    id <- which(nl$number.lineages <= n)
    dd <- setdiff(unlist(dd[id]), nl$node[id])
    
    ## select species
    ## --------------
    spec <- lapply(dd, descendants, phy = tr, type = "t", 
                   labels = TRUE, ignore.tip = TRUE)
    spec <- sapply(spec, sample, size = 1)
    y <- y[y$taxon %in% spec, ]
  }
  y
}