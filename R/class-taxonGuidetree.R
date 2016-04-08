## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2015-12-22)

## USER LEVEL CONSTRUCTOR
## ----------------------
"taxonGuidetree" <- function(ingroup, extend.ingroup = FALSE,
                             outgroup, extend.outgroup = FALSE,
                             kingdom, hybrids = FALSE,
                             tip.rank = "spec",
                             reference.rank = "auto",
                             guide.tree){
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
  
  new("taxonGuidetree", 
      ingroup = ingroup,
      extend.ingroup = extend.ingroup,
      outgroup = outgroup,
      extend.outgroup = extend.outgroup,
      kingdom = kingdom,
      hybrids = hybrids,
      tip.rank = tip.rank,
      reference.rank = reference.rank,
      guide.tree = guide.tree
  )
}

# setMethod("show",
#           signature(object = "taxon"),
#           function (object) 
#           {
#             cat("--- megaptera taxon class ---")
#             i <- object@ingroup
#             i <- paste(head(i, 2), collapse = ", ") 
#             if ( length(object@ingroup) > 2 ){
#               i <- paste(i, ", ... [", 
#                          length(object@ingroup), "]")
#             } 
#             cat("\ningroup taxon  :", i)
#             o <- object@outgroup
#             o <- paste(head(o, 2), collapse = ", ") 
#             if ( length(object@outgroup) > 2 ){
#               o <- paste(o, ", ... [", 
#                          length(object@outgroup), "]")
#             }
#             cat("\noutgroup taxon :", o)
#             cat("\nin kingdom     :", object@kingdom)
#             cat("\nhybrids        :", 
#                 ifelse(object@hybrids, "included", "excluded"))
# }
# )