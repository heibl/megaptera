## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-08-16)

setClass("taxon", 
         representation = list(
           ingroup = "list",
           extend.ingroup = "logical",
           outgroup = "list",
           extend.outgroup = "logical",
           kingdom = "character",
           hybrids = "logical",
           tip.rank = "character",
           reference.rank = "character")
)
setOldClass("phylo")
setClass(Class = "taxonGuidetree", 
         representation = list(
           guide.tree = "phylo"),
         contains = "taxon"
)

setClass("megapteraProj", 
         representation = list(
           db = "dbPars",
           taxon = "taxon",
           locus = "locus",
           align.exe = "character", 
           merge.exe = "character", 
           mask.exe = "character",
           params = "megapteraPars",
           update = "logical")
)

"megapteraProj" <- function(db, taxon, 
                            locus = locus(), 
                            align.exe = "undefined",
                            merge.exe = "undefined",
                            mask.exe = "undefined",
                            params = megapteraPars(),
                            update = FALSE){
  
  new("megapteraProj", 
      db = db,
      taxon = taxon,
      locus = locus,
      align.exe = align.exe,
      merge.exe = merge.exe,
      mask.exe = mask.exe,
      params = params,
      update = update
  )
}

setMethod("show",
          signature(object = "megapteraProj"),
          function (object) 
          {
            cat("--- megaptera project data ---")
            i <- object@taxon@ingroup
            li <- length(i)
            i <- paste(head(i, 2), collapse = ", ") 
            if ( li > 2 ){
              i <- paste(i, ", ... [", li, "]")
            } 
            cat("\ningroup taxon  :", i)
            o <- object@taxon@outgroup
            lo <- length(o)
            o <- paste(head(o, 2), collapse = ", ") 
            if ( lo > 2 ){
              o <- paste(o, ", ... [", lo, "]")
            }
            cat("\noutgroup taxon :", o)
            cat("\nin kingdom     :", object@taxon@kingdom)
            cat("\nhybrids        :", 
                ifelse(object@taxon@hybrids, "included", "excluded"))
            cat("\nlocus          :", object@locus@aliases[1])
            cat("\nexecution      :", 
                ifelse(object@params@parallel, 
                       paste("parallel on a", 
                             object@params@cluster.type,
                             "cluster with",
                             object@params@cpus, "CPUs"), 
                       "serial"))
            cat("\nupdate         :", 
                ifelse(object@update, "yes", "no"))
            cat("\nalignment      :", object@align.exe)
            cat("\nmerging        :", object@merge.exe)
            cat("\nmasking        :", object@mask.exe)
          }
)

setLocus <- function(x, locus){
#   stopifnot(inherits(x, "megapteraProj"))
#   stopifnot(inherits(locus, "locus"))
  x@locus <- locus
  x
}