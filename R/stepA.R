## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-06-12)

#' @title Step A: Creating a Project Taxonomy
#' @description Creates a project taxonomy from the NCBI taxonomy (see 
#'   \code{\link{ncbiTaxonomy}}).
#' @param x An object of class \code{\link{megapteraProj}}.
#' @return None. \code{stepA} is called for its side effects: (1) a taxonomic 
#'   classifiaction is stored in the pgSQL database; (2) a log file is written 
#'   to the \code{log} directory of the project directory.
#' @note Any subsequent call to \code{stepA} will overwrite the existing
#'   \code{taxonomy} table. As a consequence, \code{stepB} and following steps
#'   will have to be rerun also.
#' @seealso \code{\link{megapteraProj}} for bundling the project's input 
#'   information; \code{\link{stepB}} for the next step in the pipeline.
#' @export
#' @import RCurl RPostgreSQL

stepA <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  
  ## iniate logfile
  ## --------------
  logfile <- "log/stepA.log"
  if (file.exists(logfile)) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""), 
       "\nSTEP A: creating a project taxonomy from NCBI taxonomy\n",
       file = logfile)
  
  if (file.exists("ncbiTaxonomy-missing.txt")) unlink("ncbiTaxonomy-missing.txt")
  
  conn <- dbconnect(x)
  if (dbExistsTable(conn, "taxonomy"))
    slog("Existing taxonomy will be overwritten.", file = logfile)
  dbDisconnect(conn)
  
  ## get global NCBI taxonomy
  conn <- dbConnect(PostgreSQL(), dbname = "ncbitaxonomy", host = "localhost", 
                    port = 5432, user = "postgres", password = x@db@password)
  SQL <- paste("SELECT *",
               "FROM nodes",
               "JOIN names USING (id)",
               "WHERE name_class = 'scientific name'")
  tax <- dbGetQuery(conn, SQL)
  dbDisconnect(conn)
  
  ## subset to kingdom of interest
  ## -----------------------------
  tax <- taxdumpSubset(tax, mrca = x@taxon@kingdom)
  
  ## ingroup
  ## -------
  ig <- x@taxon@ingroup
  a <- sapply(ig, function(z, ncbi) any(z %in% ncbi), ncbi = tax$taxon)
  if (!all(a)){
    message("Ingroup taxa not in NCBI Taxonomy database:", paste("\n", ig[!a]))
    ig <- ig[a]
  }
  if (unique(sapply(ig, is.Linnean))){
    ingroup <- taxdumpSubset(tax, species = unlist(ig))
    ingroup <- rbind(taxdumpLineage(tax, ig[1]), ingroup) ## add lineage to root
  } else {
    ## use lapply to handle more than one taxon
    ingroup <- c(lapply(ig, taxdumpDaughters, x = tax, tip.rank = "species"),
                 lapply(ig, taxdumpLineage, x = tax))
    ingroup <- do.call(rbind, ingroup)
  }
  
  ## outgroup
  ## --------
  og <- x@taxon@outgroup
  a <- sapply(og, function(z, ncbi) any(z %in% ncbi), ncbi = tax$taxon)
  if (!all(a)){
    message("Ingroup taxa not in NCBI Taxonomy database:", paste("\n", og[!a]))
    og <- og[a]
  }
  if (unique(sapply(og, is.Linnean))){
    outgroup <- taxdumpSubset(tax, species = unlist(og))
    outgroup <- rbind(taxdumpLineage(tax, outgroup$taxon[1]), outgroup) ## add lineage to root
  } else {
    ## use lapply to handle more than one taxon
    outgroup <- c(lapply(og, taxdumpDaughters, x = tax, tip.rank = "species"),
                 lapply(og, taxdumpLineage, x = tax))
    outgroup <- do.call(rbind, outgroup)
  }
  
  tax <- unique(rbind(ingroup, outgroup))
  
  ## write parent-child taxonomy table to database
  ## ---------------------------------------------
  conn <- dbconnect(x)
  dbRemoveTable(conn, "taxonomy")
  dbWriteTable(conn, "taxonomy", tax, row.names = FALSE)
  dbSendQuery(conn, "ALTER TABLE taxonomy ADD PRIMARY KEY (id)")
  dbDisconnect(conn)
  
  slog("\n\nSTEP A finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
}
