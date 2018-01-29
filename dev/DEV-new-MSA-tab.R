

## Prepare oy dataset
library(ape)
data(woodmouse)
dna <- woodmouse
dna <- dna[order(rownames(dna)), ]
dna <- dna[1:3, 1:12]
scores <- matrix(runif(3 * 12, 0, 1), ncol = 3)
dna[2, 3:5] <- as.raw(4)
dna_ua <- del.gaps(dna)
dna_ua <- lapply(dna_ua, function(z) {dim(z) <- NULL; z})
class(dna_ua) <- "DNAbin"

## Test alignment without reliability scores
pg <- DNAbin2pg(DNAbin = dna, mode = "c")
pg2 <- pg2DNAbin(pg)
identical(dna, pg2)

## Test unaligned list of sequences
pg <- DNAbin2pg(DNAbin = dna_ua, mode = "c")
pg2 <- pg2DNAbin(pg)
identical(dna_ua, pg2)



pg <- DNAbin2pg(dna, scores)








