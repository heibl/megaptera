## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2020-01-15)

#' @title Bundle Input Data for the Pipeline
#' @description Create an object of class \code{\linkS4class{megapteraProj}} to
#'   hold all the input data necessary for running of the pipeline.
#' @param db An object of class \code{\linkS4class{dbPars}}.
#' @param taxon An object of classes \code{\linkS4class{taxon}} or
#'   \code{\linkS4class{taxonGuidetree}}.
#' @param locus An object of classes \code{\linkS4class{locus}} or
#'   \code{\linkS4class{locusRef}}.
#' @param align.exe A vector of mode \code{"character"}, giving name of the
#'   alignment program; currently only \bold{MAFFT} is allowed.
#' @param merge.exe A vector of mode \code{"character"}, giving name of the
#'   alignment merging program; currently only \bold{OPAL} is allowed.
#' @param mask.exe A vector of mode \code{"character"}, giving name of the
#'   alignment masking program; currently only \bold{Gblocks} is allowed.
#' @param params An object of class \code{\linkS4class{megapteraPars}}.
#' @param update Logical: if \code{TRUE}, the pipeline's steps are executed as
#'   if called for the first time, i.e., possibly overriding data and setting
#'   that have been previously achieved.
#' @return An object of class \code{\linkS4class{megapteraProj}}.
#' @references MAFFT: \url{http://mafft.cbrc.jp/alignment/software/}
#'
#'   OPAL: \url{http://opal.cs.arizona.edu/}
#'
#'   Gblocks: \url{http://molevol.cmima.csic.es/castresana/Gblocks.html}
#' @seealso \code{\link{dbPars}}, \code{\link{taxon}},
#'   \code{\link{taxonGuidetree}}, \code{\link{locus}}, \code{\link{locusRef}},
#'   and \code{\link{megapteraPars}} for defining of database parameters, taxa,
#'   loci, and the pipeline's parameters, respectively.
#' @include megapteraProj-class.R
#' @importFrom methods new
#' @export

"megapteraProj" <- function(db, 
                            taxon = taxon(), 
                            locus = locus(), 
                            align.exe = "undefined",
                            merge.exe = "undefined",
                            mask.exe = "undefined",
                            params = megapteraPars(),
                            update = FALSE){
  
  ## Create project file system (if it does not exist)
  ## -------------------------------------------------
  proj_path <- file.path(params@data.path, "megaptera_data/project", db@dbname)
  std_dir <- c("data", "fig", "log", "msa", "phy", "report", "scripts", "user_data")
  std_dir <- file.path(proj_path, std_dir)
  if (!dir.exists(proj_path)){
    dir.create(proj_path)
  } else {
    std_dir <- std_dir[!dir.exists(std_dir)]
  }
  if (length(std_dir)) sapply(std_dir, dir.create)
  
  ## Create instance of megapteraProj
  ## --------------------------------
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
                ifelse(object@taxon@exclude.hybrids, "excluded", "included"))
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
