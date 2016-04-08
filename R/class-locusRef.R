## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-02-20)

setOldClass("DNAbin")
setClass(Class = "locusRef", 
         representation = list(
           reference = "DNAbin"),
         contains = "locus"
)

## SET METHOD: INITIALIZE
## ----------------------
setMethod("initialize", "locusRef",
          function(.Object, kind, sql, aliases, 
                   not, adj.gene1, adj.gene2,
                   search.fields, use.genomes,
                   align.method,
                   min.identity, min.coverage,
                   reference){
            if ( !missing(aliases) ){
              .Object@kind <- kind
              .Object@sql <- sql
              .Object@aliases <- aliases
              .Object@not <- not
              .Object@adj.gene1 <- adj.gene1
              .Object@adj.gene2 <- adj.gene2
              .Object@search.fields <- search.fields
              .Object@use.genomes <- use.genomes
              .Object@align.method <- align.method
              .Object@min.identity <- min.identity
              .Object@min.coverage <- min.coverage
              .Object@reference <- reference
            }
            return(.Object)
          }
) 


## USER LEVEL CONSTRUCTOR
## ----------------------
"locusRef" <- function(..., not, 
                       search.fields = c("gene", "title"), 
                       use.genomes = TRUE,
                       align.method = "auto",
                       min.identity = 0.75,
                       min.coverage = 0.5,
                       reference, 
                       adj.gene1 = NULL, 
                       adj.gene2 = NULL,
                       check = FALSE){
  
  aliases <- c(...)
  if ( is.null(aliases) ) stop("empty '...' argument")
  
  ## determine kind of DNA
  ## to be added when 'kind' is official argument
  
  ## intergenic spacer
  ## -----------------
  if ( length(grep("IGS|intergenic spacer", aliases)) > 0 ){
    kind <- "igs"
    id <- min(grep("IGS|intergenic spacer", aliases))
    adj <- gsub(" IGS|intergenic spacer", "", aliases[id])
    adj <- unlist(strsplit(adj, "-"))
    adj.gene1 <- c(adj[1], adj.gene1)
    adj.gene2 <- c(adj[2], adj.gene2)
  } else {
    kind <- "gene" ## protein-coding gene
    if ( length(grep("[[:digit:]]S", aliases)) > 0 ) kind <- "rRNA"
    if ( length(grep("tRNA", aliases)) > 0 ) kind <- "tRNA"
    adj.gene1 <- "undefined"
    adj.gene2 <- "undefined"
  }
  
  ## GenBank uses uppercase and lowercase spelling
  ## in the same places ...
  if ( kind == "gene" ){
    aliases <- unique(c(aliases, 
                        toupper(aliases),
                        tolower(aliases), 
                        paste(toupper(substring(aliases, 1, 1)), 
                              substring(aliases, 2), sep = "")))
  }
  
  if ( missing(not) ) not <- "not"
  
  ## check if locus exists at GenBank
  ## --------------------------------
  if ( check ){
    cat("\n.. checking if locus exists on GenBank ..")
    #a <- paste("\"", aliases, "\"", sep = "")
    a <- gsub(" ", "+", aliases)
    url <- paste(a, collapse = " OR ")
    if ( !"not" %in% not ) {
      n <- paste("\"", not, "\"", sep = "")
      url <- paste(url, "NOT", not)
    }
    url <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
                 "esearch.fcgi?db=nucleotide&term=", url, 
                 "&rettype=gb&retmode=xml", sep = "")
    x <- xmlTreeParse(url, getDTD = FALSE, 
                      useInternalNodes = TRUE)
    x <- unique(xpathSApply(x, "//Count", xmlValue))
    cat("\n.. found", x[1], "records")
  }
  
  ## if GIs are given download  and
  ## assemble reference sequences
  ## ----------------------------
  if ( !inherits(reference, "DNAbin") ){
    cat("\n.. extracting reference sequence ..")
    xml <- EFetchXML(reference)
    if ( kind == "igs" ){
      reference <- lapply(reference, extractIGS, 
                          xml = xml, locus = list(adj.gene1, adj.gene2))
    } else {
      reference <- lapply(reference, extractLocus, 
                          xml = xml, locus = aliases, kind = kind)
    }
    names(reference) <- sapply(reference, names)
    failures <- sapply(reference, is.na)
    if ( all(failures) ) stop("reference sequences could not be extracted")
    if ( any (failures) ){
      reference <- reference[!failures]
      cat("\n.. success:", length(which(!failures)), "- failure:", length(which(failures)), "..")
    }
    reference <- as.DNAbin(lapply(reference, s2c))
  }
  
  new(Class = "locusRef", 
      kind = kind,
      aliases = aliases, 
      not = not, 
      sql = sql.conform(aliases[1]),
      adj.gene1 = adj.gene1,
      adj.gene2 = adj.gene2,
      search.fields = search.fields,
      use.genomes = use.genomes,
      align.method = align.method,
      min.identity = min.identity,
      min.coverage = min.coverage,
      reference = reference
  )
}

## SET METHOD: SHOW
## ----------------
setMethod("show",
          signature(object = "locusRef"),
          function (object) 
          {
            if ( object@sql == "undefined" ){
              cat("\nLocus definition: empty")
            } else {
              if ( is.matrix(object@reference) ){
                rn <- rownames(object@reference)
              } else {
                rn <- names(object@reference)
              }
              sg <- ifelse(object@use.genomes, "yes", "no")
              mi <- ifelse(object@min.identity < 0, "auto", object@min.identity)
              mc <- ifelse(object@min.coverage < 0, "auto", object@min.coverage)
              cat("\nLocus definition for", object@aliases[1],
                  "\nkind                : ", 
                  paste(object@kind, collapse = ", "),
                  "\nsearch strings      : ", 
                  paste(object@aliases, collapse = ", "),
                  "\nsearch fields       : ", 
                  paste(object@search.fields, collapse = ", "),
                  "\nuse genomes         : ", sg,
                  "\nSQL tables          : ", 
                  paste(c("acc", "spec"), 
                        object@sql, 
                        sep = "_", collapse = ", "),
                  "\nalignment method    : ", 
                  paste(object@align.method, collapse = ", "),
                  "\nminimum identity    : ", mi,
                  "\nminimum coverage    : ", mc,
                  "\nreference sequences : ", 
                  paste(rn, collapse = ", "))
            }
          }
)