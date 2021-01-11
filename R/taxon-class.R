## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-10-30)

setClass("taxon", 
         representation = list(
           ingroup = "list",
           extend.ingroup = "logical",
           outgroup = "list",
           extend.outgroup = "logical",
           kingdom = "character",
           exclude.hybrids = "logical",
           tip.rank = "character",
           reference.rank = "character")
)