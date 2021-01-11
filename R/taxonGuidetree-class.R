## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-11-13)

#' @include taxon-class.R

setOldClass("phylo")
setClass(Class = "taxonGuidetree", 
         representation = list(
           guide.tree = "phylo"),
         contains = "taxon"
)