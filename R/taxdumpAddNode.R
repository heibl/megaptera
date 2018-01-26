## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-06-12)

#' @title Utilities for NCBI Taxdump
#' @description Insert a taxon as new node to the NCBI taxonomy.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param taxon A character string giving the name of the new taxon.
#' @param rank A character string giving the rank of the new taxon.
#' @param parent A character string giving the name of the parent taxon.
#' @return None, \code{taxdumpAddNote} is called for its side effect of
#'   inserting data in the underlying PostgreSQL database.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpChildren}},
#'   \code{\link{taxdumpDropTip}}, \code{\link{taxdumpHigherRank}}, \code{\link{taxdumpMRCA}},
#'   \code{\link{taxdumpSubset}}, \code{\link{taxdump2phylo}} and
#'   \code{\link{taxdump_isTerminal}}.
#' @export

taxdumpAddNode <- function(x, taxon, rank, parent){
  
  conn <- dbconnect(x)
  
  ## check if taxon is really missing
  ## --------------------------------
  SQL <- paste("SELECT * FROM taxonomy",
               "WHERE", wrapSQL(taxon, "taxon", "="))
  SQL <- dbGetQuery(conn, SQL)
  if (nrow(SQL)) {
    message("'", taxon, "' already present in taxonomy")
    
  } else {
    ## retrieve parent data
    ## --------------------
    SQL <- paste("SELECT * FROM taxonomy",
                 "WHERE", wrapSQL(parent, "taxon", "="))
    SQL <- dbGetQuery(conn, SQL)
    if (!nrow(SQL)) stop("'", parent, "' not present in taxonomy")
    
    ## create new ID
    ## -------------
    ID <- dbGetQuery(conn, "SELECT max(id) + 1 FROM taxonomy")[1, 1]
    
    ## insert new node into database table
    ## -----------------------------------
    SQL <- paste(wrapSQL(as.integer(ID), term = NULL, boolean = NULL),
                 wrapSQL(SQL$parent_id, term = NULL, boolean = NULL),
                 wrapSQL(rank, term = NULL, boolean = NULL),
                 wrapSQL(taxon, term = NULL, boolean = NULL),
                 sep = ", ")
    SQL <- paste("INSERT INTO taxonomy", 
                 "(id, parent_id, rank, taxon)",  
                 "VALUES (", SQL, ")")
    dbSendQuery(conn, SQL)
  }
  obj <- dbDisconnect(conn)
}
