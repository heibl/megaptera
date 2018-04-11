## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-02-02)

#' @title Step B: Search and Download Sequences
#' @description For any given project (see \code{\link{megapteraProj}}), 
#' this step searches the Nucleotide database on GenBank, downloads all 
#' sequences and stores them in a postgreSQL database table.
#' @details All accessions are stored under their species name as appearing 
#' in the \emph{organism} field at GenBank, but information about infrageneric 
#' ranks is stripped off the taxon names before they are stored in the database.
#' @param x An object of class \code{\link{megapteraProj}}
#' @param update.seqs A character string determining the behaviour of 
#' \code{stepB} when it is run repeatedly on the same locus/taxon 
#' combination: \code{"no"} means that only sequences that have been made 
#' available on GenBank since the last execution \code{stepB} will be 
#' downloaded; \code{"all"} means the database table with all sequences 
#' will be removed and all sequences will be downloaded de novo. Naturely, 
#' this option has no effect if \code{stepB} is called for the first time.
#' @return None. \code{stepB} is called for its side effects: (1) strings of 
#' DNA sequences with attribute data are stored in a pgSQL database; (2) a log 
#' file is written to the current working directory.
#' @references \url{http://www.ncbi.nlm.nih.gov/books/NBK25501/}
#' @seealso \code{\link{megapteraProj}}; \code{\link{stepA}} for the preceeding 
#' and \code{\link{stepC}} for the subsequent step; \code{\link{stepBX}} for the 
#' addition of external sequences to the database.
#' @export
#' @import DBI

stepB <- function(x, update.seqs = "no"){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("undefined locus not allowed")
  if (!url.exists("https://eutils.ncbi.nlm.nih.gov"))
    stop("internet connection required for stepB")
  
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x)
  
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  dbProgress(x, "step_b", "error")
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepB.log")
  if (file.exists(logfile)) unlink(logfile)
  slog(paste("\nmegaptera", 
             packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""), 
       "\nSTEP B: searching and downloading sequences from GenBank\n",
       file = logfile)
  
  ## delete table (if it exists and is not to be updated)
  ## ----------------------------------------------------
  if (dbExistsTable(conn, acc.tab) & update.seqs == "all") {
    dbRemoveTable(conn, acc.tab)
    slog("\nRemoving TABLE", acc.tab, file = logfile)
  }
  
  ## create table (if it does not exist)
  ## -----------------------------------
  if (!dbExistsTable(conn, acc.tab)) {
    slog("\nCreating TABLE", acc.tab, file = logfile)
    SQL <- paste(acc.tab, "_pk", sep = "")
    SQL <- paste("CREATE TABLE", acc.tab, 
                 "(gi text NOT NULL,",
                 "taxon text NOT NULL,",
                 "spec_ncbi text NOT NULL,",
                 "status text,",
                 "genom text,",
                 "npos integer NOT NULL,", 
                 "identity real,",
                 "coverage real,",
                 "dna text NOT NULL,",
                 "CONSTRAINT", SQL, "PRIMARY KEY ( gi ))")
    dbSendQuery(conn, SQL)
  }
  
  ## list of species or higher taxa to be be searched for
  ## ----------------------------------------------------
  slog("\nAssembling taxon search list", file = logfile)
  ingroup <- x@taxon@ingroup
  if (unique(is.Linnean(unlist(ingroup))) & x@taxon@extend.ingroup){
    ingroup <- unique(lapply(ingroup, strip.spec))
  }
  outgroup <- x@taxon@outgroup
  if (unique(is.Linnean(unlist(outgroup))) & x@taxon@extend.outgroup){
    outgroup <- unique(lapply(outgroup, strip.spec))
  }
  search.tax <- c(ingroup, outgroup)
  
  ## search and download sequences
  ## -----------------------------
  lapply(search.tax, downloadSequences, x = x)
  
  ## declare excluded taxa: has been removed upstream to XML2acc (2016-11-03)
  
  ## rename infraspecific taxa
  ## -------------------------
  rename <- paste("SELECT taxon", 
                  "FROM", acc.tab, 
                  "WHERE status != 'excluded (indet)'")
  rename <- unique(dbGetQuery(conn, rename)$taxon)
  rename <- data.frame(gb = rename,  
                       spec = as.Linnean(rename), 
                       stringsAsFactors = FALSE)
  id <- which(rename$gb != rename$spec)
  if (length(id)){
    rename <- rename[id, ]
    rename$gb <- gsub("'", ".", rename$gb) ## Amylosporus_sp._'succulentus' # evil again!
    SQL <- paste("UPDATE", acc.tab, 
                 "SET", wrapSQL(rename$spec, term = "taxon", boolean = NULL, 
                                operator = "="),
                 "WHERE", wrapSQL(rename$gb, term = "taxon", boolean = NULL, 
                                  operator = "="))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  
  ## mark sequences too long to align in <status>
  ## --------------------------------------------
  SQL <- paste("UPDATE", acc.tab, 
               "SET status = 'excluded (too long)'",
               "WHERE npos >", x@params@max.bp)
  dbSendQuery(conn, SQL)
  
  ## Select sequences if there are > max.gi.per.spec
  ## -----------------------------------------------
  dbMaxGIPerSpec(x)
  
  ## Update relation <taxonomy>
  ## Handle species found in stepB that are not included in taxonomy table
  dbUpdateTaxonomy(x, logfile = logfile)
  
  # summary
  # -------
  total.acc <- dbGetQuery(conn, paste("SELECT count(taxon) FROM", acc.tab))
  det.acc <- dbGetQuery(conn, paste("SELECT count(taxon) FROM", acc.tab, "WHERE status !~'excluded|too'"))
  spec <- dbGetQuery(conn, paste("SELECT taxon, count(taxon) FROM", acc.tab, "WHERE status !~'excluded|too' GROUP BY (taxon)"))
  dbDisconnect(conn)
  one <- length(which(spec$count == 1))
  multiple <- length(which(spec$count > 1))
  det.spec <- one + multiple
  slog(paste("\nNumber of all accessions                              :", total.acc),
       paste("\nNumber of determined accessions                       :", det.acc),
       paste("\nNumber of determined species                          :", det.spec),
       paste("\nNumber of determined species with one accession       : ", one, 
             " (", round(one * 100 / det.spec, 1), "%)", sep = ""),
       paste("\nNumber of determined species with multiple accessions : ", multiple, 
             " (", round(multiple * 100 / det.spec, 1), "%)", sep = ""),
       file = logfile)
  
  
  slog("\n\nSTEP B finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), "\n", file = logfile)
  dbProgress(x, step = "step_b", status = "success")
}
