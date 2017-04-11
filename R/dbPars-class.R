## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-10-21)

setClass("dbPars", 
         representation = list(
           host = "character", 
           port = "numeric", 
           dbname = "character", 
           user = "character", 
           password = "character")
)