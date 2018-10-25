## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-10-25)

#' @title Utilities for NCBI Taxdump
#' @description Insert a taxon as new node to the NCBI taxonomy.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param tab A data frame in parent-child format.
#' @param rank A character string giving the rank of the new taxon, default is
#'   \code{"species"}.
#' @param taxon A character string giving the name of the new taxon.
#' @param parent A character string giving the name of the parent taxon.
#' @details \code{taxdumpAddNode} can be used in two ways: (1) If \code{x} is an
#'   object of class \code{"megapteraProj"}, the changes are made in the
#'   underlying postgreSQL database and \code{TRUE} or \code{FALSE} is returned
#'   indicating success or failure. (2) If \code{x} is a data frame in
#'   parent-child format, the changes are made to this table, which is then
#'   returned; in case of failure a warning is issued and the unchanged table is
#'   returned.
#' @return See Details section.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpChildren}},
#'   \code{\link{taxdumpDropTip}}, \code{\link{taxdumpHigherRank}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
#' @export

taxdumpAddNode <- function(x, tab, rank = "species", taxon, parent){

  #########################################
  ## MODE 1: changes in postgreSQL database
  #########################################
  if (inherits(tab, "megapteraProj")){
    
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
    
  } else {
    
    ############################################
    ## MODE 2: changes in parent-child table 'x'
    ############################################
    
    ## Assemble parent-child table if arguments are given separately
    ## -------------------------------------------------------------
    if (missing(tab)){
      if (missing(parent)) parent <- strip.spec(taxon)
      if (!parent %in% x$taxon) {
        warning("parent '", parent, "' not present in 'x'", sep = "")
        return(x)
      }
      if (rank == "species"){
        parent_rank <- x$rank[x$taxon == parent]
        # ignore "subgenus" if subgenus and genus are homonymous
        if (all(parent_rank %in% c("genus", "subgenus"))){
          parent_rank <- "genus"
        }
        if (parent_rank == "genus"){ 
          tab <- unique(data.frame(taxon = c(taxon, parent),
                                   rank = c(rank, "genus"),
                                   status = "scientific name",
                                   stringsAsFactors = FALSE))
        } else {
          tab <- unique(data.frame(taxon = c(taxon, strip.spec(taxon), parent),
                                   rank = c(rank, "genus", parent_rank),
                                   status = "scientific name",
                                   stringsAsFactors = FALSE))
        }
      } else {
        ## Taxon is not of rank "species"
        parent_rank <- x$rank[x$taxon == parent]
        tab <- data.frame(taxon = c(taxon, parent),
                         rank = c(rank, parent_rank),
                         status = "scientific name",
                         stringsAsFactors = FALSE)
      }  
    }
    
    ## Determine which separater is used by 'x' and impose it on 'tab'
    ## -----------------------------------------------------------------
    ## Beware of evil strings like 'Tuberculina sp. Ru_hy-01'
    test <- head(tab$taxon[tab$rank == "species"])
    underscore <- length(grep("[[:upper:]][[:lower:]+]_[[:lower:]]", test)) > 0
    if (underscore){
      tab$taxon <- gsub(" ", "_", tab$taxon)
    } else {
      tab$taxon <- gsub("_", " ", tab$taxon)
    }
    
    ## Check if target taxon is already present
    ## ----------------------------------------
    if (tab$taxon[1] %in% x$taxon){
      if (x@debug.level) cat("'", tab$taxon[1], "' already present in 'x'", sep = "")
      return(x)
    }
    
    ## Identify anchor taxon, i.e. the taxon of lowest rank
    ## that occurs in both 'x' and 'tab'
    ## -----------------------------------
    id <- tab$taxon %in% x$taxon
    if (!any(id)) {
      if (x@debug.level) cat("'tab' has no anchor point in 'x'")
      return(x)
    }
    ## Avoid conflict if ranks higher than anchor do not match
    id[which(id)[1]:length(id)] <- TRUE
    
    include <- which(!id)
    anchor <- tab$taxon[min(which(id))]
    
    ## Identify parent rank; necessary because 
    ## genus and subgenus can have the same names
    ## ------------------------------------------
    parent_rank <- x$rank[x$taxon == anchor]
    if (length(parent_rank) > 1){
      if (all(c("genus", "subgenus") %in% parent_rank)){
        parent_rank <- "genus"
      }
      if ("no rank" %in% parent_rank){
        parent_rank <- parent_rank[parent_rank != "no rank"]
      }
    }
  
    ## Identify anchor ID, which is the parent ID for
    ## the taxon of highest rank of the lineage that will
    ## be included
    ## --------------------------------------------------
    anchor_id <- x$id[x$taxon == anchor & x$rank == parent_rank]
    
    ## Construct new rows for taxonomy
    ## -------------------------------
    ID <- as.numeric(max(x[, c("id", "parent_id")]))
    ID <- rev(include) + ID
    tab <- data.frame(id = ID, 
                      parent_id = c(ID[-1], anchor_id), 
                      taxon = tab$taxon[include],
                      rank = tab$rank[include], 
                      status = tab$status[include],
                      stringsAsFactors = FALSE)
    
    ## Add new rows to taxonomy and return
    ## -----------------------------------
    return(rbind(x, tab))
  }
}