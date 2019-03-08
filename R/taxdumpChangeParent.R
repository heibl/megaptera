## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-03-05)

#' @title Utilities for NCBI Taxdump
#' @description Change a node's parent in the NCBI taxonomy.
#' @param x An object of class \code{\link{megapteraProj}} or a
#'   \code{data.frame} a parent-child-format.
#' @param taxon A character string giving the name of the taxon.
#' @param new.parent A character string giving the name of the new parent taxon.
#' @param orphaned.parent A character string indicating how to treat previous
#'   parent taxa, that do not have other children than \code{"taxon"}, i.e.
#'   which are orphaned.
#' @details Two choices are possible for \code{orphaned.parent}: (1)
#'   \code{"synonym"} will classify the previous, orphaned parent as a synonym
#'   of \code{"new.parent"}, whereas (2) \code{"delete"} will simply delete the
#'   orphaned parent from the taxonomy.
#' @return A \code{data.frame} in parent-child-format, if \code{class(x) ==
#'   "data.frame"}, or \code{TRUE}, if \code{class(x) == "megapteraProj"}. In
#'   the latter case \code{taxdumpChangeParent} is called for its side effect on
#'   the taxonomy table of the database.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpChildren}},
#'   \code{\link{taxdumpDropTip}}, \code{\link{taxdumpHigherRank}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
#' @importFrom DBI dbSendQuery   
#' @export

taxdumpChangeParent <- function(x, taxon, new.parent, orphaned.parent = "synonym"){
  
  if (inherits(x, "megapteraProj")){
    tax <- dbReadTaxonomy(x)
  } else {
    tax <- x
  }
  orphaned.parent <- match.arg(orphaned.parent, c("delete", "synonym"))
  
  if (!taxon %in% tax$taxon) {
    stop("Taxon '", taxon, "' not present in taxonomy", sep = "")
  }
  
  ## Look up IDs
  ## -----------
  opid <- tax$parent_id[tax$taxon == taxon] ## old parent's ID
  oname <- tax$taxon[tax$id == opid] ## old parent's name
  npid <- tax$id[tax$taxon == new.parent] ## new parent's ID
  sid <- tax$id[tax$parent_id == opid & tax$taxon != taxon] ## sister taxa's ID

  if (inherits(x, "megapteraProj")){
    conn <- dbconnect(x)
    SQL <- paste("UPDATE taxonomy",
                 "SET", wrapSQL(npid, "parent_id", "="), 
                 "WHERE", wrapSQL(taxon, "taxon", "="))
    dbSendQuery(conn, SQL)
    cat("Parent of '", taxon, "' changed from '", tax$taxon[tax$id == opid],
        "' [", opid, "] to '", new.parent, "' [", npid, "]", sep = "")
    
    
    if (!length(sid)){
      if (orphaned.parent == "delete"){
        SQL <- paste("DELETE FROM taxonomy",
                     "WHERE", wrapSQL(opid, "id", "="))
        dbSendQuery(conn, SQL)
        cat("\nPrevious parent ('", oname, "' [", opid, "]) deleted", sep = "")
      } else {
        ## orphaned.parent == "synonym"
        SQL <- paste("UPDATE taxonomy",
                     "SET status='synonym',", 
                     wrapSQL(npid, "id", "="),
                     "WHERE", wrapSQL(opid, "id", "="))
        dbSendQuery(conn, SQL)
        cat("\nStatus of previous parent ('", oname, "' [", opid, " -> ", npid, "])",
            " changed to 'synonym'", sep = "")
      }
    }
    dbDisconnect(conn)
    cat("\nNo output - changes affect only the PostgreSQL database")
    invisible(TRUE)
    
  } else {
    
    ## change and return data frame 'tax'
    ## ----------------------------------
    tax$parent_id[tax$taxon == taxon] <- npid
    cat("Parent of '", taxon, "' changed from '", tax$taxon[tax$id == opid],
        "' [", opid, "] to '", new.parent, "' [", npid, "]", sep = "")
    if (!length(sid)){
      if (orphaned.parent == "delete"){
        tax <- tax[tax$id != opid, ]
        cat("\nPrevious parent ('", oname, "' [", opid, "]) deleted", sep = "")
      } else {
        ## orphaned.parent == "synonym"
        tax$status[tax$id == opid] <- "synonym"
        tax$id[tax$id == opid] <- npid
        cat("\nStatus of previous parent ('", oname, "' [", opid, " -> ", npid, "])",
            " changed to 'synonym'", sep = "")
      }
    }
    cat("\nChanges affect  only the returned data frame (and not the PostgreSQL database)")
    return(tax)
  } 
}
