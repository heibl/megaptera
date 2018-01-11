## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-01-11)

#' @title Class "megapteraProj"
#' @description S4 class for holding all input data of the megaptera pipeline.
#' @slot db An object of class \code{\linkS4class{dbPars}}.
#' @slot taxon An object of classes \code{\linkS4class{taxon}} or
#'   \code{\linkS4class{taxonGuidetree}}.
#' @slot locus An object of classes \code{\linkS4class{locus}} or
#'   \code{\linkS4class{locusRef}}.
#' @slot align.exe A vector of mode \code{"character"}, giving name of the
#'   alignment program; currently only \bold{MAFFT} is allowed.
#' @slot merge.exe A vector of mode \code{"character"}, giving name of the
#'   alignment merging program; currently only \bold{OPAL} is allowed.
#' @slot mask.exe A vector of mode \code{"character"}, giving name of the
#'   alignment masking program; currently only \bold{Gblocks} is allowed.
#' @slot params An object of class \code{\linkS4class{megapteraPars}}.
#' @slot update Logical: if \code{TRUE}, the pipeline's steps are executed as
#'   if called for the first time, i.e., possibly overriding data and setting
#'   that have been previously achieved.
#' @seealso \code{\link{dbPars}}, \code{\link{taxon}}, \code{\link{locus}}, and
#'   \code{\link{megapteraPars}} for defining of database parameters, taxa,
#'   loci, and the pipeline's parameters, respectively.
#' @include dbPars.R dbPars-class.R taxon-class.R taxonGuidetree-class.R

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




