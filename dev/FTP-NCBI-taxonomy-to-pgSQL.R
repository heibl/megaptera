##############################################################
# TASK        : download NCBI taxdump database from FTP server
#               and create ncbitaxonomy database in PostgreSQL
# AUTHOR      : Christoph Heibl
# LAST CHANGE : 2017-02-16
##############################################################
library(RPostgreSQL)
setwd("/Users/heibl/Documents/r/phylogeny")

## download and decompress taxdump
## -------------------------------
system("curl -OL ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", 
       ignore.stderr = TRUE, wait = TRUE)
untar("taxdump.tar.gz", exdir = "taxdump")

conn <- dbConnect(PostgreSQL(), dbname = "ncbitaxonomy", host = "localhost", port = 5432,
                  user = "postgres", password = "oxalis")

## 1. nodes: the hierarchical structure of the classification
## ----------------------------------------------------------
nodes <- read.table("taxdump/nodes.dmp", sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)
names(nodes) <- c("id", "parent_id", "rank", "embl_code", "division_id",
                  "inherited_div_flag", "genetic_code_id", "inherited_gc_flag",
                  "mitochondrial_genetic_code_id", "inherited_MGC_flag", "GenBank_hidden_flag",
                  "hidden_subtree_root_flag", "comments")
nodes <- nodes[, 1:3]

dbRemoveTable(conn, "nodes")
SQL <- paste("CREATE TABLE public.nodes(",
             "id integer NOT NULL,",
             "parent_id integer,",
             "rank text,",
             "CONSTRAINT nodes_pkey PRIMARY KEY (id))")
dbSendQuery(conn, SQL)
dbWriteTable(conn, "nodes", nodes, append = TRUE, row.names = FALSE)
remove(nodes)

## 2. names: taxon names
## ---------------------
# taxnames <- read.table("taxdump/names.dmp", sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)
taxnames <- scan("taxdump/names.dmp", sep = "\n", what = "c", strip.white = TRUE)
taxnames <- lapply(taxnames, function(z) unlist(strsplit(z, "\t[|]\t|\t[|]")))
id <- sapply(taxnames, length)
if (all(id == 4)) taxnames <- do.call(rbind, taxnames)
taxnames <- as.data.frame(taxnames, stringsAsFactors = FALSE)
names(taxnames) <- c("id", "taxon", "unique_name", "name_class")

dbRemoveTable(conn, "names")
SQL <- paste("CREATE TABLE public.names(",
  "id integer NOT NULL,",
  "taxon text NOT NULL,",
  "unique_name text NOT NULL,",
  "name_class text NOT NULL,",
  "CONSTRAINT names_pkey PRIMARY KEY (id, taxon, unique_name, name_class))")
dbSendQuery(conn, SQL)
dbWriteTable(conn, "names", taxnames, append = TRUE, row.names = FALSE)
remove(taxnames)

dbDisconnect(conn)


