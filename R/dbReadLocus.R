## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-01-04)

dbReadLocus <- function(x, provenance = ".", tag){
  
  tip.rank <- x@taxon@tip.rank
  conn <- dbconnect(x)

  ## select columns: curently not used
  ## ---------------------------------
  #   col.set <- c("all", "gb", "sel")
  #   cols <- match.arg(cols, col.set, several.ok = TRUE)
  #   if ( "all" %in% cols ) cols <- col.set[-1]
  #   cols  <- paste("_", cols, "$", sep = "")
  #   cols <- paste(cols, collapse = "|")
  #   loc <- loc[, grep(cols, colnames(loc))]
  #   
  
  ## names of acc_* tables
  ## ---------------------
  acc <- dbTableNames(x, "acc")
#   not.empty <- paste("SELECT count(gi) FROM", acc)
#   not.empty <- sapply(not.empty, dbGetQuery, conn = conn)
#   not.empty <- unlist(not.empty) > 0
#   acc <- acc[not.empty]
  
  ## names of spec_* or gen_* tables
  ## -------------------------------
  msa <- dbTableNames(x, tip.rank)
  
  ## gb:
  dbGenBank <- function(conn, tab, tip.rank, provenance, tag){
    tag <- ifelse(missing(tag), "",
                  paste("WHERE", wrapSQL(tag, "t.tag")))
    SQL <- paste("SELECT taxon",
                 "FROM", tab,
                 "WHERE", wrapSQL(provenance, term = "genom"))
    SQL <- paste("SELECT t.xxx, count(a.taxon)",
                 "AS", gsub("acc", "gb", tab),
                 "FROM taxonomy AS t", 
                 "LEFT JOIN (", SQL, ") AS a", 
                 "ON t.spec = a.taxon",
                 tag,
                 "GROUP BY t.xxx",
                 "ORDER BY t.xxx")
    SQL <- gsub("xxx", tip.rank, SQL)
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
    if ( tip.rank == "gen" ) tag <- ""
    SQL <- ifelse(tip.rank == "spec", 
                  "FROM taxonomy AS t",
                  "FROM (SELECT DISTINCT gen FROM taxonomy) AS t")
    SQL <- paste("SELECT t.xxx, s.status", 
                 "AS", gsub(tip.rank, "sel", tab),
                  SQL,
                 "LEFT JOIN", tab, "AS s",
                 "ON t.xxx = s.xxx", 
                 tag,
                 "ORDER BY t.xxx")
    SQL <- gsub("xxx", tip.rank, SQL)
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
  if ( !is.null(tab.sel) ){
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