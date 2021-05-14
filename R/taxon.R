## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-11-24)

#' @include taxon-class.R
#' @importFrom methods new
#' @export

## USER LEVEL CONSTRUCTOR
## ----------------------
"taxon" <- function(ingroup, extend.ingroup = FALSE,
                    outgroup, extend.outgroup = FALSE,
                    kingdom, exclude.hybrids = FALSE,
                    tip.rank = "species",
                    reference.rank = "auto"){
  
  if (missing(ingroup)){
    new("taxon", 
        ingroup = list("undefined"),
        extend.ingroup = extend.ingroup,
        outgroup = list("undefined"),
        extend.outgroup = extend.outgroup,
        kingdom = "undefined",
        exclude.hybrids = exclude.hybrids,
        tip.rank = tip.rank,
        reference.rank = reference.rank)
  } else {
    
    ##
    if (!taxdumpSanity(ingroup, quiet = TRUE)){
      stop("inconsistent 'ingroup' object: run taxdumpSanity(ingroup) for diagnosis")
    }
    if (!taxdumpSanity(outgroup, quiet = TRUE)){
      stop("inconsistent 'outgroup' object: run taxdumpSanity(outgroup) for diagnosis")
    }
   
    tip.rank <- match.arg(tip.rank, c("genus", "species"))
    
    new("taxon", 
        ingroup = ingroup,
        extend.ingroup = extend.ingroup,
        outgroup = outgroup,
        extend.outgroup = extend.outgroup,
        kingdom = kingdom,
        exclude.hybrids = exclude.hybrids,
        tip.rank = tip.rank,
        reference.rank = reference.rank)
  }
}

setMethod("show",
          signature(object = "taxon"),
          function (object) {
            
            arg <- c("tip rank", "ingroup", "is extended" ,  
                     "outgroup", "is extended",    
                     "in kingdom", "exclude.hybrids", "guide tree")
            arg <- format(arg)
            
            formatTaxa <- function(taxa){
              root <- taxa$taxon[taxa$id == taxa$parent_id]
              spec <- taxa$taxon[taxa$rank == "species" & taxa$status == "scientific name"]
              n <- length(spec)
              spec <- sort(spec)[c(1, n)]
              spec <- paste(spec, collapse = ifelse(n > 2, ", ..., ", ", "))
              paste0(n, " species of ", root, ": ", spec)
            }
            out <- c(
              object@tip.rank,
              formatTaxa(object@ingroup), 
              ifelse(object@extend.ingroup, "yes", "no"),
              formatTaxa(object@outgroup),
              ifelse(object@extend.outgroup, "yes", "no"),
              object@kingdom,
              ifelse(object@exclude.hybrids, "excluded", "included"),
              ifelse(inherits(object, "taxonGuidetree"), 
                     "user-defined", "taxonomy-based")
            )
            out <- paste(arg, out, sep = " : ")
            out <- c("--- megaptera taxon class ---", out)
            out <- paste("\n", out, sep = "")
            cat(paste(out, collapse = ""))
          }
)