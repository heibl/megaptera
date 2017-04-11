## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-11-03)

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