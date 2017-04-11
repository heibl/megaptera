# TASK: retreive a taxon's classification from ncbitaxonomy (pgSQL)
# AUTHOR: Christoph Heibl
# LAST UPDATE: 2017-01-10

library(megaptera)
library(RPostgreSQL)

conn <- dbConnect(PostgreSQL(), dbname = "ncbitaxonomy", host = "localhost", 
                  port = 5432, user = "postgres", password = "oxalis")
SQL <- paste("SELECT *",
             "FROM nodes",
             "JOIN names USING (id)",
             "WHERE name_class = 'scientific name'")
system.time(tt <- dbGetQuery(conn, SQL))
dbDisconnect(conn)

## 1: Simple tests:
system.time(p <- taxdump2phylo(tt, taxon = "Viperidae"))
system.time(p <- taxdump2phylo(tt, taxon = "Cantharellales", tip.rank = "species"))
system.time(p <- taxdump2phylo(tt, taxon = "Basidiomycota", tip.rank = "species"))

# plot(p, show.tip.label = FALSE, no.margin = TRUE, cex = 0.5)
# plot(ladderize(p), show.tip.label = FALSE, no.margin = TRUE, cex = 0.5)

## 2: Test with orders of Basidiomycota
b <- read.tree("/Users/heibl/Documents/r/tol/results/basidiomycota-guidetree.tre")
bo <- b$tip.label
p <- lapply(bo, taxdump2phylo, x = tt)
names(p) <- bo


## graft subtrees onto guide tree
## ------------------------------
for ( i in bo ){
  if (inherits(p[[i]], "phylo")){
    b <- bind.tree(b, p[[i]], which(b$tip.label == i))
  } else {
    b$tip.label[b$tip.label == i] <- p[[i]]
  }
}

plot(b, show.tip.label = FALSE, no.margin = TRUE, cex = 0.5)

## save
## ----
write.tree(b, "/Users/heibl/Documents/r/phylogeny/basidiomycota/data/guidetree.tre")
b <- read.tree("/Users/heibl/Documents/r/phylogeny/basidiomycota/data/guidetree.tre")

