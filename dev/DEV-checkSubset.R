library(megaptera)
library(biodiv)
project.name <- "coleoptera"
load(paste0("data/", project.name, ".megapteraProj"))

subset <- refMatch(ref = "coleoptera", format = "megaptera")



z <- checkSubset(x, subset)

zsm <- read.fas("data/ColeopteraZSM20180122.fas")
zsm2 <- read.FASTA("data/ColeopteraZSM20180122.fas")
zsm <- names(zsm2)
zsm <- strsplit(zsm, split = "|", fixed = TRUE)
zsm <- sapply(zsm, function(z) z[2])
zsm <- unique(zsm)

zsm <- read.fas("data/test.fas")
miss <- z$taxon[z$present == "no"]

## from NCBI and ZSM
miss <- setdiff(miss, zsm)

