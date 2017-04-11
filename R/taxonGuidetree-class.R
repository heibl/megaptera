## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-11-03)

#' @include taxon-class.R

setOldClass("phylo")
setClass(Class = "taxonGuidetree", 
         representation = list(
           guide.tree = "phylo"),
         contains = "taxon"
)