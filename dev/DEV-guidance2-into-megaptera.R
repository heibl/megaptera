## Integrate GUIDANCE2 into megaptera
library(megaptera)
library(polenta)
library(stringr)
library(parallel)
library(doSNOW)

project.name <- "chiroptera-europe"
load(paste0("data/", project.name, ".megapteraProj"))
load("data/mt.rda")
x <- setLocus(x, mt[["cox1"]])
gene <- x@locus@sql
msa.tab <- paste(x@taxon@tip.rank, gsub("^_", "", gene), sep = "_")
seqs <- dbReadDNA(x, msa.tab, regex = TRUE, blocks = "ignore")
## arguments
sequences <- del.gaps(seqs)
bootstrap <- 3
msa.exec <- "/usr/local/bin/mafft"
ncore <- 24
method <- "auto"
res1 <- guidance2(sequences, bootstrap = 3, 
                  msa.exec = "/usr/local/bin/mafft", ncore = 24)

dim(seqs@scores)

res2 <- guidance(sequences, bootstrap = 3, 
                  msa.exec = "/usr/local/bin/mafft", ncore = 24)
dim(res2@scores)
