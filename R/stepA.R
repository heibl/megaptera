## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-11-20)

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
#'   information and \code{\link{ncbiTaxonomy}} for downloading the NCBI
#'   taxonomy; the next step in the pipeline is \code{\link{stepB}}.
#' @export
#' @import RCurl RPostgreSQL
#' @importFrom utils write.table

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
    slog("Existing taxonomy will be overwritten", file = logfile)
  dbDisconnect(conn)
  
  ## Get global NCBI taxonomy
  ## ------------------------
  conn <- dbConnect(PostgreSQL(), dbname = "ncbitaxonomy", host = "localhost", 
                    port = 5432, user = "postgres", password = x@db@password)
  SQL <- paste("SELECT *",
               "FROM nodes",
               "JOIN names USING (id)",
               "WHERE name_class = 'scientific name'")
  tax <- dbGetQuery(conn, SQL)
  dbDisconnect(conn)
  
  ## Subset to kingdom of interest
  ## -----------------------------
  tax <- taxdumpSubset(tax, mrca = x@taxon@kingdom)
  
  ## INGROUP
  ## -------
  ig <- x@taxon@ingroup
  ## Check if ingroup is present in NCBI taxonomy
  ## --------------------------------------------
  a <- sapply(ig, function(z, ncbi) any(z %in% ncbi), ncbi = tax$taxon)
  if (x@params@debug.level){
    if (!all(a)){
      missing_ig <- sort(unlist(ig[!a]))
      slog("\n", length(missing_ig), " ingroup taxa not in NCBI Taxonomy database:", 
           paste("\n", missing_ig), sep = "", file = logfile)
      if (x@params@debug.level > 1){
        write.table(missing_ig, file = "data/ingroup-not-NCBI-taxonomy.txt",
                    quote = FALSE, row.names = FALSE, col.names = "INGROUP" )
        slog("\nThe list of missing ingroup taxa was written to", 
             "'data/ingroup-not-NCBI-taxonomy.txt'", file = logfile)
      }
    } else {
      slog("\nAll ingroup taxa are present in NCBI Taxonomy database", file = logfile)
    }
  }
  ig <- ig[a]
  
  if (unique(sapply(unlist(ig), is.Linnean))){
    ingroup <- taxdumpSubset(tax, species = unlist(ig))
    ingroup <- rbind(taxdumpLineage(tax, ig[1]), ingroup) ## add lineage to root
  } else {
    ## use lapply to handle more than one taxon
    ingroup <- c(lapply(ig, taxdumpDaughters, tax = tax, tip.rank = "species"),
                 lapply(ig, taxdumpLineage, tax = tax))
    ingroup <- do.call(rbind, ingroup)
  }
  
  ## OUTGROUP
  ## --------
  og <- x@taxon@outgroup
  ## Check if ingroup is present in NCBI taxonomy
  ## --------------------------------------------
  a <- sapply(ig, function(z, ncbi) any(z %in% ncbi), ncbi = tax$taxon)
  if (x@params@debug.level){
    if (!all(a)){
      missing_og <- sort(unlist(og[!a]))
      slog("\n", length(missing_og), " outgroup taxa not in NCBI Taxonomy database:", 
           paste("\n", og[!a]), sep = "", file = logfile)
      if (x@params@debug.level > 1){
        write.table(missing_og, file = "data/outgroup-not-NCBI-taxonomy.txt",
                    quote = FALSE, row.names = FALSE, col.names = "OUTGROUP" )
        slog("\nThe list of missing outgroup taxa was written to", 
             "'data/outgroup-not-NCBI-taxonomy.txt'", file = logfile)
      }
    } else {
      slog("\nAll outgroup taxa are present in NCBI Taxonomy database", file = logfile)
    }
  }
  og <- og[a]
  if (unique(sapply(unlist(og), is.Linnean))){
    outgroup <- taxdumpSubset(tax, species = unlist(og))
    outgroup <- rbind(taxdumpLineage(tax, outgroup$taxon[1]), outgroup) ## add lineage to root
  } else {
    ## use lapply to handle more than one taxon
    outgroup <- c(lapply(og, taxdumpDaughters, tax = tax, tip.rank = "species"),
                 lapply(og, taxdumpLineage, tax = tax))
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
