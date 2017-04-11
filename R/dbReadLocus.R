## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-02-01)

#' @importFrom DBI dbDisconnect dbGetQuery
#' @export

dbReadLocus <- function(megProj, provenance = ".", tag){
  
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
    SQL <- ifelse(tip.rank == "genus", "regexp_replace(taxon, ' .+$', '') AS taxon", "taxon")
    SQL <- paste("SELECT", SQL, "FROM", tab)
    SQL <- paste("SELECT t.taxon, count(a.taxon)",
                 "AS", gsub("acc", "gb", tab),
                 "FROM taxonomy AS t", 
                 "LEFT JOIN (", SQL, ") AS a", 
                 "ON t.taxon = a.taxon",
                 "WHERE", wrapSQL(tip.rank, "t.rank", "="),
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
  dbSelected <- function(conn, tab, tip.rank, tag){
    tag <- ifelse(missing(tag), "",
                  paste("WHERE", wrapSQL(tag, "t.tag")))
    if (tip.rank == "genus") tag <- ""
    SQL <- paste("(SELECT taxon",
                  "FROM taxonomy",
                  "WHERE", wrapSQL(tip.rank, "rank", "="))
    SQL <- paste("SELECT t.taxon, s.status", 
                 "AS", gsub(tip.rank, "sel", tab),
                 "FROM", SQL,
                 ") AS t LEFT JOIN", tab, "AS s",
                 "ON t.taxon = s.taxon", 
                 tag,
                 "ORDER BY t.taxon")
    tab <- dbGetQuery(conn, SQL)
    tab[is.na(tab)] <- 0
    rownames(tab) <- tab[, 1]
    tab[, 2, drop = FALSE]
  }
  tab.sel <- lapply(msa, dbSelected, conn = conn, 
                    # tag = tag, 
                    tip.rank = tip.rank)
  tab.sel <- do.call(cbind, tab.sel)
  
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
  obj
}