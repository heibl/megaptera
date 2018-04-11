## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-02-21)

#' @title Utilities for NCBI Taxdump
#' @description Insert a taxon as new node to the NCBI taxonomy.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param tab A data frame in parent-child format.
#' @param rank A character string giving the rank of the new taxon, default is
#'   \code{"species"}.
#' @param taxon A character string giving the name of the new taxon.
#' @param parent A character string giving the name of the parent taxon.
#' @return Logical, \code{TRUE} if the taxon was already present or was added
#'   successfully into the NCBI taxonomy, \code{FALSE} if the taxon could not be
#'   inserted due to a missing anchor point. Obviously, \code{taxdumpAddNote} is
#'   called for its side effect of inserting data in the underlying PostgreSQL
#'   database.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpChildren}},
#'   \code{\link{taxdumpDropTip}}, \code{\link{taxdumpHigherRank}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
#' @export

taxdumpAddNode <- function(x, tab, rank = "species", taxon, parent){
  
  if (!missing(tab)){
    id <- which(tab$rank == rank)
    taxon <- tab$taxon[id]
    parent <- tab$taxon[(id + 1):(nrow(tab))]
    names(parent) <- tab$rank[(id + 1):(nrow(tab))]
  }
  
  conn <- dbconnect(x)
  
  ## Check if taxon is really missing
  ## --------------------------------
  SQL <- paste("SELECT * FROM taxonomy",
               "WHERE", wrapSQL(taxon, "taxon", "="))
  SQL <- dbGetQuery(conn, SQL)
  if (nrow(SQL)) {
    dbDisconnect(conn)
    message("'", taxon, "' already present in taxonomy")
    return(TRUE)
  }
  
  ## Retrieve parent data. This can be done over a vector of
  ## ranks that are ordered from lower to higher (i.e. towards the root)
  ## The parent of lowest available rank is then used as an anchor point
  ## -------------------------------------------------------------------
  SQL <- paste("SELECT * FROM taxonomy",
               "WHERE", wrapSQL(parent, "taxon", "=", NULL))
  SQL <- lapply(SQL, dbGetQuery, conn = conn)
  id <- sapply(SQL, nrow)
  if (sum(id) == 0) {
    message("parent taxa not present in taxonomy")
    return(FALSE)
  }
  
  ## Construct new lineage to insert
  ## -------------------------------
  ID <- "SELECT max(parent_id), max(id) FROM taxonomy"
  ID <- dbGetQuery(conn, ID)
  ID <- max(ID)
  anchor <- which(id == 1)[1]
  if (anchor == 1){
    ## Case 1: taxon can be inserted directly, e.g. genus is already present
    ## ---------------------------------------------------------------------
    obj <- data.frame(id = ID + 1, parent_id = SQL[[anchor]]$id, 
                      rank, taxon, stringsAsFactors = FALSE)
  } else {
    ## Case 2: taxon cannot be inserted directly, e.g. genus is absent
    ## ---------------------------------------------------------------
    ID <- rev((1:anchor) + ID)
    obj <- data.frame(id = ID, 
                      parent_id = c(ID[-1], SQL[[anchor]]$id), 
                      rank = tab$rank[1:anchor], 
                      taxon = tab$taxon[1:anchor], 
                      stringsAsFactors = FALSE)
  }
  
  ## Insert new node into database table
  ## -----------------------------------
  SQL <- paste(wrapSQL(obj$id, term = NULL, boolean = NULL),
               wrapSQL(obj$parent_id, term = NULL, boolean = NULL),
               wrapSQL(obj$rank, term = NULL, boolean = NULL),
               wrapSQL(obj$taxon, term = NULL, boolean = NULL),
               sep = ", ")
  SQL <- paste("INSERT INTO taxonomy", 
               "(id, parent_id, rank, taxon)",  
               "VALUES (", SQL, ")")
  lapply(SQL, dbSendQuery, conn = conn)
  
  
  dbDisconnect(conn)
  return(TRUE)
}
