## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-04-20)

#' @title NCBI Taxonomy Database
#' @description Retrieve taxonomic classification from the taxonomy database at
#'   the National Center for Biotechnology Information (NCBI).
#' @param x Database connection parameters, either as object of class 
#'   \code{\link{dbPars}} or \code{\link{megapteraProj}}.
#' @details The NCBI taxonomy database will be downloaded via FTP in "taxdump" 
#'   format, unpacked, and translated into a data frame object. In a second 
#'   step, the data frame table will be stored in a postgreSQL database called 
#'   \code{"ncbitaxonomy"}.
#' @return \code{ncbiTaxonomy} is called for its side effect (see Details).
#' @references NCBI Taxonomy Database website:
#'   \url{http://www.ncbi.nlm.nih.gov/taxonomy}
#' @references Federhen, Scott. 2012. The NCBI taxonomy database. \emph{Nucleic
#'   Acids Research} \bold{40}: DI36-DI43.
#' @seealso \code{\link{stepA}} for extracting the taxonomic information 
#'   relevant for a given \bold{megaptera} project.
#'   \code{\link{dbUpdateTaxonomy}} and \code{\link{dbReadTaxonomy}} for storing
#'   and retrieving taxonomic information in a \bold{megaptera} project
#'   database.
#' @importFrom DBI dbDisconnect dbRemoveTable dbSendQuery dbWriteTable
#' @export

ncbiTaxonomy <- function(x){
  
  ## check input object
  ## ------------------
  if (!inherits(x, c("dbPars", "megapteraProj"))){
    stop("x is not of classes 'dbPars' or 'megapteraProj'")
  }
  if (inherits(x, "megapteraProj")){
    x <- x@db
  }
  
  ## check if database exists ...
  ## ----------------------------
  conn <- dbConnect(RPostgreSQL::PostgreSQL(),
                    host = x@host,
                    port = x@port,
                    user = x@user,
                    password = x@password)
  sql <- paste("SELECT 1 FROM pg_database WHERE",
               sql.wrap("ncbitaxonomy", term = "datname"))
  if ( nrow(dbGetQuery(conn, sql)) == 1 ){
    cat("\ndatabase 'ncbitaxonomy' exists ans will be updated")  
  } else {
    ## .. and create if it does not exist
    ## ----------------------------------
    cat("\ndatabase 'ncbitaxonomy' created") 
    sql <- paste("CREATE DATABASE ncbitaxonomy",
                 "WITH ENCODING='UTF8'",
                 "CONNECTION LIMIT=-1;")
    dbSendQuery(conn, sql)
  }
  dbDisconnect(conn)
  
  ## connect to 'ncbitaxonomy'
  ## -------------------------
  conn <- dbConnect(RPostgreSQL::PostgreSQL(),
                    host = x@host,
                    port = x@port,
                    user = x@user,
                    password = x@password,
                    dbname = "ncbitaxonomy")
  
  ## download and decompress taxdump
  ## -------------------------------
  unlink("taxdump.tar.gz")
  unlink("taxdump", recursive = TRUE)
  system("curl -OL ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", 
         ignore.stderr = TRUE, wait = TRUE)
  untar("taxdump.tar.gz", exdir = "taxdump")
  
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
  
  ## disconnect and eliminate traces
  ## -------------------------------
  dbDisconnect(conn)
  unlink("taxdump.tar.gz")
  unlink("taxdump", recursive = TRUE)
}