## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2021-03-16)

#' @title Utilities for NCBI Taxdump
#' @description Change a node's parent in the NCBI taxonomy.
#' @param x An object of class \code{\link{megapteraProj}} or a
#'   \code{data.frame} a parent-child-format.
#' @param old.taxon A character string giving the name of the taxon that will
#'   be replaced.
#' @param new.taxon A character string giving the new name of for the taxon in
#'   question.
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
#' @importFrom crayon %+% bold cyan magenta silver
#' @importFrom DBI dbSendQuery   
#' @export

taxdumpRename <- function(x, old.taxon, new.taxon){
  
  if (inherits(x, "megapteraProj")){
    ## Need to set drop.extinct = FALSE, in order not to loose potential parentials
    ## that just have been added
    tax <- dbReadTaxonomy(x, drop.extinct = FALSE)
  } else {
    tax <- x
  }
  
  if (!old.taxon %in% tax$taxon) {
    stop("old taxon name '", old.taxon, "' not present in taxonomy", sep = "")
  }
  if (new.taxon %in% tax$taxon) {
    stop("new taxon name '", new.taxon, "' already in use", sep = "")
  }
  
  id <- tax$id[tax$taxon == old.taxon]
  
  ## Make changes
  if (inherits(x, "megapteraProj")){
    
    conn <- dbconnect(x)
    SQL <- paste("UPDATE taxonomy",
                 "SET", wrapSQL(new.taxon, "taxon", "="), 
                 "WHERE", wrapSQL(old.taxon, "taxon", "="))
    dbSendQuery(conn, SQL)
    cat(silver("Taxon name '" %+% magenta$bold(old.taxon) %+% "' [" %+% cyan(id) %+% 
                "] changed to '" %+% magenta$bold(new.taxon) %+% "\n", sep = ""))
    cat(silver("No output - changes affect only the PostgreSQL database\n"))
    invisible(TRUE)
    
  } else {
    
    ## change and return data frame 'tax'
    ## ----------------------------------
    tax$taxon[tax$taxon == old.taxon] <- new.taxon
    cat(silver("Taxon name '" %+% magenta$bold(old.taxon) %+% "' [" %+% cyan(id) %+% 
                 "] changed to '" %+% magenta$bold(new.taxon) %+% "\n", sep = ""))
    cat(silver("Changes affect only the returned data frame (and not the PostgreSQL database)"))
    return(tax)
  } 
}
