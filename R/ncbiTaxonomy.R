## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-10-09)

#' @title NCBI Taxonomy Database
#' @description Retrieve taxonomic classification from the Taxonomy database
#'   maintained by the National Center for Biotechnology Information (NCBI).
#' @param x Database connection parameters, either as object of class
#'   \code{\link{dbPars}} or \code{\link{megapteraProj}}. Alternatively, can be
#'   a character string giving the path to the megaptera file system (see
#'   \code{\link{megapteraInit}}).
#' @param quiet Logical, indicating if diagnostic messages should be printed on
#'   screen.
#' @details The NCBI taxonomy database will be downloaded via FTP in "taxdump"
#'   format, unpacked, and translated into data frames. In a second step, the data
#'   frames will be stored in a postgreSQL database called
#'   \code{"ncbitaxonomy"}. Any existing data will be overwritten in the
#'   process.
#'
#'   There is no versioning of the Taxonomy database like the release numbers of
#'   GenBank. The taxonomy FTP dump files are updated hourly so it might be
#'   reasonable to update them together with GenBank sequences data.
#' @return \code{ncbiTaxonomy} is called for its side effect (see Details).
#' @references NCBI Taxonomy Database website:
#'   \url{http://www.ncbi.nlm.nih.gov/taxonomy}
#' @references Federhen, Scott. 2012. The NCBI taxonomy database. \emph{Nucleic
#'   Acids Research} \bold{40}: DI36-DI43.
#' @seealso \code{\link{stepA}} for extracting the taxonomic information
#'   relevant for a given megaptera project; \code{\link{dbUpdateTaxonomy}} and
#'   \code{\link{dbReadTaxonomy}} for storing and retrieving taxonomic
#'   information in a megaptera project database;
#'   \code{\link{dbSummaryTaxonomy}} for a short numerical summary of the
#'   project taxonomy.
#' @importFrom DBI dbConnect dbDisconnect dbRemoveTable dbSendQuery dbWriteTable
#' @importFrom tools md5sum
#' @importFrom utils download.file read.table untar
#' @export

ncbiTaxonomy <- function(x, quiet = FALSE){
  
  ## Check input object
  ## ------------------
  if (inherits(x, c("dbPars", "megapteraProj"))){
    if (inherits(x, "megapteraProj")){
      path <- file.path(x@params@data.path, "megaptera_data/gb_taxonomy")
      x <- x@db
    }
  } else {
    path <- file.path(x, "megaptera_data/gb_taxonomy")
    x <- NULL
    if (!dir.exists(path)){
      stop("x must be either the path to the megaptera file system ", 
           "or of classes 'dbPars' or 'megapteraProj'")
    }
  }
  
  ## 1. Determine if existing version of taxdump is still up to date
  ## ---------------------------------------------------------------
  new_version <- TRUE ## default
  unlink(file.path(path, "taxdump.tar.gz.md5"))
  download.file("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5",
                file.path(path, "taxdump.tar.gz.md5"))
  md5 <- scan(file.path(path, "taxdump.tar.gz.md5"), what = "c", quiet = TRUE)[1]
  if (file.exists(file.path(path, "INFO"))){
    info_old <- scan(file.path(path, "INFO"), sep = "\n", what = "c", quiet = TRUE)
    info_old <- sapply(info_old, strsplit, split = ": ")
    info_old <- do.call(rbind, info_old)
    md5_old <- info_old[info_old[, 1] == "MD5", 2]
    if (!quiet) message("\nFound 'taxdump' download from ", 
                        info_old[info_old[, 1] == "Timestamp", 2])
    new_version <- md5_old != md5
    if (new_version){
      if (!quiet) message("\nThere is a new 'taxdump' version available")
    } else {
      if (!quiet) message("\n'taxdump' version is still up to date ... ")
    }
  } 
  
  ## 2. Download taxdump via FTP and check file with MD5
  ## ---------------------------------------------------
  if (new_version){
    
    if (!quiet) message("\nDownloading 'taxdump' via FTP ... ")
    unlink(file.path(path, "taxdump.tar.gz"))
    unlink(file.path(path, "taxdump"), recursive = TRUE)
    
    download.file("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz",
                  file.path(path, "taxdump.tar.gz"))
    
    md5 <- scan(file.path(path, "taxdump.tar.gz.md5"), what = "c", quiet = TRUE)[1]
    if (md5sum(file.path(path, "taxdump.tar.gz")) != md5){
      stop("MD5 checksums do not match")
    }
    
    info <- c(paste("Timestamp:", as.character(Sys.time())),
              paste("MD5:", md5),
              paste("MegapteraVersion:", packageDescription("megaptera")$Version))
    info <- paste(info, collapse = "\n")
    write(info, file.path(path, "INFO"))
  }
  if (!quiet) message("done")
  
  if (is.null(x)) return() ## stop here if only path was given
  
  ## 3. Decompress taxdump
  ## ---------------------
  if (!quiet) message("\nDecompressing 'taxdump.tar.gz' ... ")
  untar(file.path(path, "taxdump.tar.gz"), exdir = file.path(path, "taxdump"))
  if (!quiet) message("done")
  
  ## Check if database exists ...
  ## ----------------------------
  conn <- dbConnect(RPostgreSQL::PostgreSQL(),
                    host = x@host,
                    port = x@port,
                    user = x@user,
                    password = x@password)
  sql <- paste("SELECT 1 FROM pg_database WHERE",
               wrapSQL("ncbitaxonomy", "datname", "="))
  if (nrow(dbGetQuery(conn, sql)) == 1){
    if (!quiet) cat("\nDatabase 'ncbitaxonomy' exists and will be updated")  
  } else {
    ## .. and create if it does not exist
    ## ----------------------------------
    if (!quiet) cat("\nDatabase 'ncbitaxonomy' created") 
    sql <- paste("CREATE DATABASE ncbitaxonomy",
                 "WITH ENCODING='UTF8'",
                 "CONNECTION LIMIT=-1;")
    dbSendQuery(conn, sql)
  }
  dbDisconnect(conn)
  
  ## Connect to 'ncbitaxonomy'
  ## -------------------------
  conn <- DBI::dbConnect(RPostgreSQL::PostgreSQL(),
                         host = x@host,
                         port = x@port,
                         user = x@user,
                         password = x@password,
                         dbname = "ncbitaxonomy")
  
  ## 1. nodes: the hierarchical structure of the classification
  ## ----------------------------------------------------------
  if (!quiet) cat("\nReading table 'nodes' ... ")
  nodes <- read.table(file.path(path, "taxdump/nodes.dmp"), sep = "|", 
                      strip.white = TRUE, stringsAsFactors = FALSE)
  names(nodes) <- c("id", "parent_id", "rank", "embl_code", "division_id",
                    "inherited_div_flag", "genetic_code_id", "inherited_gc_flag",
                    "mitochondrial_genetic_code_id", "inherited_MGC_flag", "GenBank_hidden_flag",
                    "hidden_subtree_root_flag", "comments")
  nodes <- nodes[, 1:3]
  
  if (!quiet) cat("done\nWriting table 'nodes' to postgreSQL database 'ncbitaxonomy' ... ")
  dbRemoveTable(conn, "nodes")
  SQL <- paste("CREATE TABLE public.nodes(",
               "id integer NOT NULL,",
               "parent_id integer,",
               "rank text,",
               "CONSTRAINT nodes_pkey PRIMARY KEY (id))")
  dbSendQuery(conn, SQL)
  dbWriteTable(conn, "nodes", nodes, append = TRUE, row.names = FALSE)
  remove(nodes)
  if (!quiet) cat("done")
  
  ## 2. names: taxon names
  ## ---------------------
  if (!quiet) cat("\nReading table 'names' ... ")
  # taxnames <- read.table("taxdump/names.dmp", sep = "|", strip.white = TRUE, stringsAsFactors = FALSE)
  taxnames <- scan(file.path(path, "taxdump/names.dmp"), sep = "\n", 
                   what = "c", strip.white = TRUE, quiet = TRUE)
  taxnames <- lapply(taxnames, function(z) unlist(strsplit(z, "\t[|]\t|\t[|]")))
  id <- sapply(taxnames, length)
  if (all(id == 4)) taxnames <- do.call(rbind, taxnames)
  taxnames <- as.data.frame(taxnames, stringsAsFactors = FALSE)
  names(taxnames) <- c("id", "taxon", "unique_name", "name_class")
  if (!quiet) cat("done")
  
  ## Check for duplicate entries
  d <- duplicated(taxnames)
  if (any(d)){
    dd <- taxnames[d, ]
    if (!quiet) cat("\nRemoving", nrow(dd), "duplicate", ifelse(nrow(dd) == 1, "entry","entries"),
                    "from 'names':", formatDF(dd))
    taxnames <- taxnames[!d, ]
  }
  
  if (!quiet) cat("\nWriting table 'names' to postgreSQL database 'ncbitaxonomy' ... ")
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
  if (!quiet) cat("done")
  
  ## Disconnect and eliminate traces
  ## -------------------------------
  if (!quiet) cat("\nCleaning up ... ")
  dbDisconnect(conn)
  # unlink(file.path(path, "taxdump.tar.gz"))
  unlink(file.path(path, "taxdump.tar.gz.md5"))
  unlink(file.path(path, "taxdump"), recursive = TRUE)
  # if (!quiet) cat("done")
}