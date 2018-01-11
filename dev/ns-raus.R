library(ggplot2)
library(ips)

dna <- read.fas("/Users/heibl/Documents/r/decay_evo_8_11_2017/msa_1218x37466_14-08-2017.fas")
dna2 <- deleteEmptyCells(dna)
iupac <- c(n = 240, "-" = 4)
n <- as.raw(iupac)
nn <- apply(dna, 2, function(z) length(which(z %in% n)))
nn <- nn/nrow(dna)

qs <- read.table("/Users/heibl/Documents/r/decay_evo_8_11_2017/guidance_score-37466_14-08-2017.txt")
qs <- unlist(qs[1, , drop = TRUE])


plot(qs ~ nn)

hist(nn)
hist(qs)

## IUPAC ambiguity code
## --------------------

