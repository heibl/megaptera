## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-12-06)

#' @include locus-class.R
#' @importFrom methods new
#' @export
  
## USER LEVEL CONSTRUCTOR
## ----------------------
"locus" <- function(..., not,
                    search.fields = c("gene", "title"), 
                    use.genomes = TRUE,
                    align.method = "auto",
                    min.identity = 0.75,
                    min.coverage = 0.5,
                    check = FALSE){
  if ( missing(...) ){
    new(Class = "locus", 
        kind = "undefined",
        aliases = "undefined", 
        not = "undefined", 
        adj.gene1 = "undefined",
        adj.gene2 = "undefined",
        sql = "undefined",
        search.fields = search.fields,
        use.genomes = use.genomes,
        align.method = align.method,
        min.identity = min.identity,
        min.coverage = min.coverage
    )
  } else {
    aliases <- c(...)
    
    ## intergenic spacer
    ## -----------------
    if ( length(grep("intergenic spacer", aliases)) > 0 ){
      kind <- "igs"
      adj <- gsub(" intergenic spacer", "", aliases[1])
      adj <- unlist(strsplit(adj, "-"))
      adj.gene1 <- adj[1]
      adj.gene2 <- adj[2]
    } else {
      kind <- ifelse(length(grep("[[:digit:]]S", aliases)) > 0, 
                     "rRNA", "gene")
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
      cat("\nChecking if locus exists on GenBank ..")
      #a <- paste("\"", aliases, "\"", sep = "")
      a <- gsub(" ", "+", aliases)
      url <- paste(a, collapse = " OR ")
      if ( !"not" %in% not ) {
        n <- paste("\"", not, "\"", sep = "")
        url <- paste(url, "NOT", not)
      }
      x <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
                  "esearch.fcgi?db=nucleotide&term=", url, 
                  "&rettype=gb&retmode=xml")
      x <- xmlTreeParse(getURL(x), getDTD = FALSE, 
                        useInternalNodes = TRUE)
      x <- unique(xpathSApply(x, "//Count", xmlValue))
      cat("\n.. found", x[1], "records")
    }
    new(Class = "locus", 
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
        min.coverage = min.coverage
    )
  }
}

## SET METHOD: SHOW
## ----------------
setMethod("show",
          signature(object = "locus"),
          function (object) 
          {
            if ( object@sql == "undefined" ){
              cat("\nLocus definition: empty")
            } else {
              mi <- ifelse(object@min.identity < 0, "auto", object@min.identity)
              mc <- ifelse(object@min.coverage < 0, "auto", object@min.coverage)
              cat("\nLocus definition for", object@aliases[1],
                  "\nkind                : ", 
                  paste(object@kind, collapse = ", "),
                  "\nsearch strings      : ", 
                  paste(object@aliases, collapse = ", "),
                  "\nsearch fields       : ", 
                  paste(object@search.fields, collapse = ", "),
                  "\nuse genomes         : ", 
                  paste(object@use.genomes, collapse = ", "),
                  "\nSQL tables          : ", 
                  paste(c("acc", "spec"), 
                        object@sql, 
                        sep = "_", collapse = ", "),
                  "\nalignment method    : ", 
                  paste(object@align.method, collapse = ", "),
                  "\nminimum identity    : ", mi,
                  "\nminimum coverage    : ", mc)
            }
          }
)