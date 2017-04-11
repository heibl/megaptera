## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-11-03)

#' @include taxon-class.R taxonGuidetree-class.R

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




