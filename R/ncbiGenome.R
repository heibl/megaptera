## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2018-04-11)

## to do:
## extend outgroup to genus if it is a single species

#'@title Find Organelle Genomes on NCBI GenBank
#'@description Finds the most representative organelle genomes of a given taxon
#'  for use as reference sequences.
#'@param x An object of class \code{\link{megapteraProj}}.
#'@param organelle A character string, either \code{"mitochondrion"},
#'  \code{"chloroplast"}, or any unambiguous abreviation of these.
#'@param mrca A vector of mode \code{"character"}, can be \code{"ingroup"},
#'  \code{"outgroup"}, or both. In the latter case, reference genomes are
#'  searched for in all taxa descending from the MRCA of ingroup and outgroup.
#'@param n Numeric, the maximum number of genomes that will be chosen. Depending
#'  on the classification of the taxon as returned by \code{\link{stepA}}, the
#'  actual number of genomes returned can be less than \code{n}.
#'@details \code{ncbiGenome} uses a four-step algorithm to produce a
#'  taxonomically balanced sample of reference organelle genomes: \enumerate{
#'  \item Determine the root taxon for both ingroup and outgroup. \item Find all
#'  organelle genomes present on NCBI GenBank for this taxa. \item Using the
#'  taxonomic classifiaction, find the \code{n - x} basal lineages of the entire
#'  set of genomes; thereby \code{x} is often greater than 0 depending on the
#'  branching pattern (topology) encoded by the classifiaction. \item For each
#'  lineage, randomly choose one organelle genome and return the results as data
#'  frame (see \code{Value} section). }
#'@return a data frame with three columns: \item{taxon}{scientific name as Latin
#'  binomial} \item{gb}{UID: GenBank number} \item{gi}{alternative UID (GIs will
#'  no longer be supported after august 2016!)}
#'@references NCBI Orgenelle Genome Resources:
#'  \url{http://www.ncbi.nlm.nih.gov/genome/organelle/}
#'@seealso \code{\link{locusRef}} to set reference sequences.
#'@export
#'@importFrom RCurl url.exists
#'@import XML

ncbiGenome <- function(x, organelle, mrca = c("ingroup", "outgroup"), n = 5){
  
  ## CHECKS
  ## ------
  if (!url.exists("https://eutils.ncbi.nlm.nih.gov"))
    stop("internet coection required for ncbiGenome")
  
  ## set organelle
  ## -------------
  fn <- ifelse(dir.exists("results"), "log/ncbiGenome.log", "")
  organelle <- match.arg(organelle, c("mitochondrion", "chloroplast"))
  slog("\n.. organelle     :", organelle, file = fn)
  ## tag can either be 'mitochondrion' or 'mitochondrial DNA'
  organelle <- gsub("drion", "dri*", organelle)
  
  ## find common root of ingroup and outgroup
  ## ----------------------------------------
  if (all(c("ingroup", "outgroup") %in% mrca)) mrca <- "both"
  common.root <- findRoot(x, mrca)
  # common.root <- head(common.root$taxon, 1) ## slugs (2017-10-11)
  
  core <- function(this.root, organelle, fn){
    ## query ENTREZ and save history on server
    ## ---------------------------------------
    xml <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/", 
                  "eutils/esearch.fcgi?",
                  "tool=megaptera",
                  "&email=heibsta@gmx.net",
                  "&usehistory=y",
                  "&retmax=9999",
                  "&db=nucleotide",
                  "&term=(", this.root, 
                  "[Organism]+AND+", organelle,
                  "[Title]+AND+complete+genome[Title])"
    )
    # cat(xml)
    
    ## get and parse results via eFetch from history server
    ## ----------------------------------------------------
    xml <- robustXMLparse(xml, logfile = fn)
    webEnv <- xpathSApply(xml, fun = xmlToList,
                          path = "//eSearchResult/WebEnv")
    queryKey <- xpathSApply(xml, fun = xmlToList,
                            path = "//eSearchResult/QueryKey")
    nn <- as.numeric(xpathSApply(xml, fun = xmlValue,
                                 path = "//eSearchResult/Count"))
    list(webEnv = webEnv, queryKey = queryKey, nn = nn)
  }
  
  ## search iteratively beginning from the MRCA down to the root-of-life
  ## -------------------------------------------------------------------
  for (i in 1:nrow(common.root)){
    
    this.root <- common.root$taxon[i] ## slugs (2017-10-11)
    slog("\n.. common root   :", this.root, file = fn) 
    out <- core(this.root, organelle, fn)
    if (!out$nn){
      slog(" - no genomes available")
    } else {
      break
    }
  }
  
  ## loop over sliding window
  ## ------------------------
  retmax <- 50
  sw <- seq(from = 0, to = out$nn, by = retmax)
  sw <- data.frame(from = sw, to = c(sw[-1] - 1, out$nn))
  slog("\n.. posting", out$nn, "UIDs on Entrez History Server ..", file = fn)
  b <- ifelse(nrow(sw) == 1, "batch", "batches")
  slog("\n.. retrieving full records in", nrow(sw), b, "..\n", file = fn)
  
  y <- data.frame()
  for (i in 1:nrow(sw)) {
    
    ## get XML with full records
    ## -------------------------
    xml <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
                  "esummary.fcgi?tool=megaptera&email=heibsta@gmx.net",
                  "&db=nucleotide&query_key=", out$queryKey, 
                  "&WebEnv=", out$webEnv,
                  "&retmode=xml",
                  "&retstart=", sw$from[i], 
                  "&retmax=", retmax)
    
    ## parse XML: this step is error-prone and therefore
    ## embedded into try()
    ## -------------------
    xml <- robustXMLparse(xml, logfile = fn)
    # saveXML(xml, "aaa.xml"); system("open -t aaa.xml")
    if (is.null(xml)) next
    
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
  if (length(indet)){
    y <- y[-indet, ]
  }
  y$taxon <- strip.infraspec(y$taxon)
  y$taxon <- gsub(" ", "_", y$taxon)
  y <- y[!duplicated(y$taxon), ]
  
  ## create a phylogenetically balanced sample
  ## of genomes in case the number of available
  ## genomes exceeds n (the number of desired genomes)
  ## -------------------------------------------------
  if (nrow(y) > n){
    
    y <- y[sample(1:nrow(y), n), ]
    # tr <- ncbiTaxonomy(as.list(y$taxon), 
    #                    kingdom = x@taxon@kingdom,
    #                    megapteraProj = x)
    # tr <- tax2tree(tr)
    # tr <- compute.brlen(tr)
    # # tr <- cophenetic.phylo(tr)
    # branching.order <- names(sort(branching.times(tr), decreasing = TRUE))
    # branching.order <- as.numeric(branching.order)
    # dd <- lapply(branching.order, descendants, phy = tr, type = "d")
    # nl <- sapply(dd, length)
    # nl[-1] <- nl[-1] - 1
    # nl <- data.frame(node = branching.order,
    #                  number.lineages = cumsum(nl))
    # id <- which(nl$number.lineages <= n)
    # dd <- setdiff(unlist(dd[id]), nl$node[id])
    # 
    # ## select species
    # ## --------------
    # spec <- lapply(dd, descendants, phy = tr, type = "t", 
    #                labels = TRUE, ignore.tip = TRUE)
    # spec <- sapply(spec, sample, size = 1)
    # y <- y[y$taxon %in% spec, ]
  }
  y
}