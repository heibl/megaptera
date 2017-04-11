## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-12-06)

#' @include megapteraProj-class.R
#' @importFrom methods new
#' @export

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
