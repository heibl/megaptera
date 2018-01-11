## develop taxon nodelabels

tax <- dbReadTaxonomy(x, tip.rank = "species", root = "mrca")
tax <- taxdumpSubset(tax, "Ciidae", root = "mrca")
tr <- taxdump2phylo(tax, tip.rank = "species")
tr <- ladderize(tr)
tr <- fixNodes(tr)
plot(tr)



plot(tr, no.margin = TRUE)
edgelabels(obj, id, cex = 0.5, adj = c(0.5, 1.1), frame = "n")
