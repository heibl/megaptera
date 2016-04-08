## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-02-02)

## SET METHOD: INITIALIZE: see class-megapteraProj.R

## USER LEVEL CONSTRUCTOR
## ----------------------
"taxon" <- function(ingroup, extend.ingroup = FALSE,
                    outgroup, extend.outgroup = FALSE,
                    kingdom, hybrids = FALSE,
                    tip.rank = "spec",
                    reference.rank = "auto"){
 
  ingroup <- unique(ingroup); outgroup <- unique(outgroup)
  if ( is.factor(ingroup) ) ingroup <- levels(ingroup)[ingroup]
  if ( is.factor(outgroup) ) outgroup <- levels(outgroup)[outgroup]
  if ( is.character(ingroup) ) ingroup <- as.list(ingroup)
  if ( is.character(outgroup) ) outgroup <- as.list(outgroup)
  
  ## enforce canonical reference rank name
  ## -------------------------------------
  if ( reference.rank != "auto" ){
    reference.rank <- sqlTaxonomyHeader(reference.rank)
  }
  
  new("taxon", 
      ingroup = ingroup,
      extend.ingroup = extend.ingroup,
      outgroup = outgroup,
      extend.outgroup = extend.outgroup,
      kingdom = kingdom,
      hybrids = hybrids,
      tip.rank = tip.rank,
      reference.rank = reference.rank
  )
}

setMethod("show",
          signature(object = "taxon"),
          function (object) {
            
            arg <- c("ingroup", "is extended" ,  
                     "outgroup", "is extended",    
                     "in kingdom", "hybrids", "guide tree")
            arg <- format(arg)
            
            formatTaxa <- function(taxa){
              n <- length(taxa)
              taxa <- paste(head(taxa, 2), collapse = ", ") 
              if ( n > 2 ) taxa <- paste(taxa, ", ... [", n, "]")
              taxa
            }
            out <- c(formatTaxa(object@ingroup), 
                     ifelse(object@extend.ingroup, "yes", "no"),
                     formatTaxa(object@outgroup),
                     ifelse(object@extend.outgroup, "yes", "no"),
                     object@kingdom,
                     ifelse(object@hybrids, "included", "excluded"),
                     ifelse(inherits(object, "taxonGuidetree"), 
                            "user-defined", "taxonomy-based")
                     )
            out <- paste(arg, out, sep = " : ")
            out <- c("--- megaptera taxon class ---", out)
            out <- paste("\n", out, sep = "")
            cat(paste(out, collapse = ""))
}
)