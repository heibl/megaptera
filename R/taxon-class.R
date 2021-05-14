## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2021-05-11)

setClass("taxon", 
         representation = list(
           ingroup = "data.frame",
           extend.ingroup = "logical",
           outgroup = "data.frame",
           extend.outgroup = "logical",
           kingdom = "character",
           exclude.hybrids = "logical",
           tip.rank = "character",
           reference.rank = "character")
)