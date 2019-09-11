## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-04-15)

#' @title Taxon-Locus-Crosstable
#' @description Create a data frame that contains for each taxon and locus the
#'   number of accessions found on GenBank/BOLD/etc. and if the taxon was
#'   selected for alignment.
#' @param megProj An object of class
#'   \code{\linkS4class{megapteraProj}}.
#' @param provenance A characters string as used in the column \code{provenance}
#'   the \code{acc_<locus>} tables in the postgreSQL database. Depending on your
#'   actual pipeline setup these may be \code{"ncbi"}, \code{"ncbi annotated"},
#'   or \code{"bold"} as well as any user-defined strings.
#' @param tag \emph{Deprecated!} A characters string as used in the column
#'   \code{tag} in the \code{taxonomy} table in the postgreSQL database.
#'   Depending on your choices for parameter values in \code{\link{taxon}} or
#'   \code{\link{taxonGuidetree}} these may be \code{"ingroup (NCBI)"},
#'   \code{"extended ingroup (NCBI)"}, \code{"outgroup (NCBI)"}, or
#'   \code{"extended outgroup (NCBI)"} as well as any user-defined strings.
#' @param subset A vector of mode \code{"character"}, that can be used to choose
#'   a subset of the total taxa available. 
#' @return A data frame with two types of colums
#' \item{gb_<locus>}{gives the number of accessions per taxon and locus}
#' \item{sel_<locus>}{states if a taxon was selected for final alignment}
#' @seealso \code{\link{checkSpecLocus}}
#' @importFrom DBI dbDisconnect dbGetQuery
#' @export

dbReadLocus <- function(megProj, provenance = ".", tag, subset){
  
  tip.rank <- megProj@taxon@tip.rank
  conn <- dbconnect(megProj)
  
  ## names of acc_* tables
  ## ---------------------
  acc <- dbTableNames(megProj, "acc")

  ## names of spec_* or gen_* tables
  ## -------------------------------
  msa <- dbTableNames(megProj, tip.rank)
  
  ## gb:
  dbGenBank <- function(conn, tab, tip.rank, provenance, tag){
    # tag <- ifelse(missing(tag), "",
    #               paste("AND", wrapSQL(tag, "t.tag")))
    SQL <- ifelse(tip.rank == "genus", 
                  "regexp_replace(taxon, ' .+$', '') AS taxon", 
                  "taxon")
    SQL <- paste("SELECT", SQL, "FROM", tab)
    SQL <- paste("SELECT t.taxon, count(a.taxon)",
                 "AS", gsub("acc", "gb", tab),
                 "FROM taxonomy AS t", 
                 "LEFT JOIN (", SQL, ") AS a", 
                 "ON t.taxon = a.taxon",
                 "WHERE", wrapSQL(tip.rank, "t.rank", "="),
                 "AND t.status='scientific name'",
                 # tag,
                 "GROUP BY t.taxon",
                 "ORDER BY t.taxon")
    tab <- dbGetQuery(conn, SQL)
    rownames(tab) <- tab[, 1]
    tab[, 2, drop = FALSE]
  }
  tab.gb <- lapply(acc, dbGenBank, conn = conn, tip.rank = tip.rank,
                   # tag = tag,
                   provenance = provenance)
  tab.gb <- do.call(cbind, tab.gb)
  
  ## selected sequences:
  ## -------------------
  selected_loci <- "SELECT DISTINCT locus FROM species_sequence ORDER BY locus"
  selected_loci <- dbGetQuery(conn, selected_loci)$locus
  if (!is.null(selected_loci)){
    out <- gsub("__", "_", paste("sel", selected_loci, sep = "_"))
    SQL <- paste("(SELECT taxon, status",
                 "FROM species_sequence",
                 "WHERE", wrapSQL(selected_loci, "locus", "=", NULL), ")")
    SQL <- paste("SELECT t.taxon, s.status", 
                 "AS", out,
                 "FROM", SQL, "AS s",
                 "FULL JOIN taxonomy AS t ON t.taxon = s.taxon", 
                 "WHERE", wrapSQL(tip.rank, "t.rank", "="),
                 "AND t.status='scientific name'",
                 "ORDER BY t.taxon")
    
    
    tab.sel <- lapply(SQL, dbGetQuery, conn = conn)
    taxon <- tab.sel[[1]][, "taxon"]
    tab.sel <- lapply(tab.sel, function(z) z[, 2, drop = FALSE])
    tab.sel <- do.call(cbind, tab.sel)
    tab.sel[is.na(tab.sel)] <- 0
    rownames(tab.sel) <- taxon
  } else {
    tab.sel <- NULL
  }
  
  
  dbDisconnect(conn)
  
  ## create complete table
  ## ---------------------
  if (!is.null(tab.sel)){
    obj <- cbind(tab.gb, tab.sel)
    ## order columns
    loc <- rep(gsub("acc_", "", acc), each = 2)
    loc <- paste(c("gb", "sel"), loc, sep = "_")
    loc <- loc[loc %in% names(obj)]
    obj <- obj[, loc]
  } else {
    obj <- tab.gb
  }
  
  ## Create subset (optional)
  ## ------------------------
  if (!missing(subset)){
    subset <- gsub("_", " ", subset)
    subset <- intersect(subset, rownames(obj))
    if (!length(subset)) {
      stop("subset is emtpy")
    }
    obj <- obj[subset, ]
  }
  obj
}