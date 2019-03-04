## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-03-04)

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
  SQL <- paste("SELECT parent_id, id, taxon, rank, name_class AS status",
               "FROM nodes",
               "JOIN names USING (id)",
               "WHERE name_class = 'scientific name'",
               "OR name_class = 'synonym'"
               )
  tax <- dbGetQuery(conn, SQL)
  dbDisconnect(conn)
  
  ## Subset to kingdom of interest
  ## -----------------------------
  tax <- taxdumpSubset(tax, mrca = x@taxon@kingdom)
  
  ## INGROUP
  ## -------
  ig <- x@taxon@ingroup
  
  ## This is rahter dirty, but should be a quick fix for the discrepancy that
  ## arises between ICZN and Fauna Europea, which did not adopt its Gender Agreement
  ## -------------------------------------------------------------------------------
  renderGender <- function(b){
    fem <- grep("a$", b)
    mas <- grep("um$", b)
    b <- c(b, gsub("a$", "um", b[fem]), gsub("a$", "us", b[fem]),
           gsub("um$", "a", b[mas]))
    hyphen <- grep("-", b) ## NCBI: Satyrium walbum (instead of Satyrium w-album)
    unique(c(b, gsub("-", "", b[hyphen])))
  }
  ig <- lapply(ig, renderGender)
  
  ## Check if/which ingroup taxa are present in NCBI taxonomy
  ## --------------------------------------------------------
  a <- sapply(ig, function(z, ncbi) any(z %in% ncbi), ncbi = tax$taxon)
  if (x@params@debug.level){
    if (!all(a)){
      missing_ig <- ig[!a]
      slog("\n", length(missing_ig), " ingroup taxa not in NCBI Taxonomy database:", 
           formatSpecList(missing_ig), sep = "", file = logfile)
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
  
  ## Adjust accepted names/synonyms as user-defined
  ## ----------------------------------------------
  # tax2 <- tax
  ## ig[which(sapply(ig, function(x, y) x %in% y, x = "Hyloicus pinastri"))]
  for (i in seq_along(ig)[]){
    cat("\n", i, " ")
    tax <- taxdumpSynonym(tax, binomials = ig[[i]], keep.acc = FALSE, 
                   quiet = FALSE, keep.syn = TRUE)
    # if (!taxdumpSanity(tax)) break
  }
  
  if (unique(sapply(unlist(ig), is.Linnean))){
    ingroup <- taxdumpSubset(tax, species = sapply(ig, head, n = 1))
    ingroup <- rbind(taxdumpLineage(tax, ig[1]), ingroup) ## add lineage to root
  } else {
    ## Use lapply to handle more than one taxon
    ingroup <- c(lapply(ig, taxdumpChildren, tax = tax, tip.rank = "species"),
                 lapply(ig, taxdumpLineage, tax = tax))
    ingroup <- do.call(rbind, ingroup)
  }
  
  ## OUTGROUP
  ## --------
  og <- x@taxon@outgroup
  ## Check if ingroup is present in NCBI taxonomy
  ## --------------------------------------------
  a <- sapply(og, function(z, ncbi) any(z %in% ncbi), ncbi = tax$taxon)
  if (x@params@debug.level){
    if (!all(a)){
      missing_og <- og[!a]
      slog("\n", length(missing_og), " outgroup taxa not in NCBI Taxonomy database:", 
           formatSpecList(missing_og), sep = "", file = logfile)
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
  
  ## Adjust accepted names/synonyms as user-defined
  ## ----------------------------------------------
  for (i in seq_along(og)){
    cat("\n", i, " ")
    tax <- taxdumpSynonym(tax, binomials = og[[i]], quiet = FALSE, 
                          keep.acc = FALSE, keep.syn = FALSE)
  }
  
  if (unique(sapply(unlist(og), is.Linnean))){
    outgroup <- taxdumpSubset(tax, species = sapply(og, head, n = 1))
    outgroup <- rbind(taxdumpLineage(tax, outgroup$taxon[1]), outgroup) ## add lineage to root
  } else {
    ## use lapply to handle more than one taxon
    outgroup <- c(lapply(og, taxdumpChildren, tax = tax, tip.rank = "species"),
                 lapply(og, taxdumpLineage, tax = tax))
    outgroup <- do.call(rbind, outgroup)
  }
  
  tax <- unique(rbind(ingroup, outgroup))
  
  ## Check if taxonomy is sane
  ## -------------------------
  test <- taxdumpSanity(tax)
  if (!test){
    stop("debug me!")
  }
  
  ## Write parent-child taxonomy table to database
  ## ---------------------------------------------
  ## Note: 'id' has to be part of the PK because synonyms can be homonyms
  ## (by different authors, e.g. Calothysanis). Hence, the check for 
  ## duplicated accepted names has to be done upstream by tamxdumpSanity()
  ## ---------------------------------------------
  conn <- dbconnect(x)
  dbRemoveTable(conn, "taxonomy")
  dbWriteTable(conn, "taxonomy", tax, row.names = FALSE)
  dbSendQuery(conn, "ALTER TABLE taxonomy ADD PRIMARY KEY (id, taxon, rank)")
  dbDisconnect(conn)
  
  slog("\n\nSTEP A finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
}
