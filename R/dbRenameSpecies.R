## This code is part of the megaptera package
## © C. Heibl 2015 (last update: 2018-01-22)

#' @title Rename a Species
#' @description Rename a species in all tables of the database.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param from A vector of mode \code{"character"} giving the species name that should be replaced.
#' @param to A vector of mode \code{"character"} giving the new species name.
#' @return None, \code{dbRenameSpecies} is called for its side effect.
#' @seealso \code{\link{dbDeleteSpec}}, \code{\link{dbExcludeSpec}}, \code{\link{dbIncludeSpec}}
#' @importFrom DBI dbDisconnect dbSendQuery
#' @export

dbRenameSpecies <- function(x, from, to){
 
  conn <- dbconnect(x)
  
  ## Step 1: Taxonomy table
  ## ----------------------
  tax <- paste("UPDATE taxonomy",
               "SET", wrapSQL(to, "taxon", operator = "="),
               "WHERE", wrapSQL(from, "taxon", operator = "="))
  dbSendQuery(conn, tax)

  
  ## Step 2: Accession tables
  ## ------------------------
  acc <- dbTableNames(x, tabs = "acc")
  if (length(acc)){
    SQL <- paste("UPDATE", acc, 
                 "SET", wrapSQL(to, "taxon", operator = "="),
                 "WHERE", wrapSQL(from, "taxon", operator = "="))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  
  ## Step 3: MSA tables
  ## ------------------------
  msa <- dbTableNames(x, tabs = x@taxon@tip.rank)
  if (length(msa)){
    SQL <- paste("UPDATE", msa, 
                 "SET", wrapSQL(to, "taxon", operator = "="),
                 "WHERE", wrapSQL(from, "taxon", operator = "="))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  
  dbDisconnect(conn)
}