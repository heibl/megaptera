## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-11-03)

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

