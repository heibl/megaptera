## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2021-05-11)

#' @include taxonGuidetree-class.R
#' @importFrom methods new
#' @export

## USER LEVEL CONSTRUCTOR
## ----------------------
"taxonGuidetree" <- function(ingroup, extend.ingroup = FALSE,
                             outgroup, extend.outgroup = FALSE,
                             kingdom, exclude.hybrids = FALSE,
                             tip.rank = "species",
                             reference.rank = "auto",
                             guide.tree){
  if (!taxdumpSanity(ingroup, quiet = TRUE)){
    stop("inconsistent 'ingroup' object: run taxdumpSanity(ingroup) for diagnosis")
  }
  if (!taxdumpSanity(outgroup, quiet = TRUE)){
    stop("inconsistent 'outgroup' object: run taxdumpSanity(outgroup) for diagnosis")
  }
  tip.rank <- match.arg(tip.rank, c("genus", "species"))
  
  new("taxonGuidetree", 
      ingroup = ingroup,
      extend.ingroup = extend.ingroup,
      outgroup = outgroup,
      extend.outgroup = extend.outgroup,
      kingdom = kingdom,
      exclude.hybrids = exclude.hybrids,
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
#             cat("\nexclude.hybrids        :", 
#                 ifelse(object@exclude.hybrids, "included", "excluded"))
# }
# )