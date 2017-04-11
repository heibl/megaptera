## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-02-16)

setClass("megapteraPars", 
         representation = list(
           parallel = "logical",
           cpus = "numeric",
           cluster.type = "character",
           update.seqs = "character", # step B
           retmax = "numeric", # step B
           max.gi.per.spec = "numeric", # step B
           max.bp = "numeric", # step B
           reference.max.dist = "numeric", # step C
           min.seqs.reference = "numeric", # step D
           fract.miss = "numeric", # step G
           filter1 = "numeric",  # step G
           filter2 = "numeric", # step G
           filter3 = "numeric", # step G
           filter4 = "numeric",  # step G
           block.max.dist = "numeric", # step G
           min.n.seq = "numeric", # step G
           max.mad = "numeric", # step H
           gb1 = "numeric", # step F + G
           gb2 = "numeric", # step F + G
           gb3 = "numeric", # step F + G
           gb4 = "numeric", # step F + G
           gb5 = "character") # step F + G
)
