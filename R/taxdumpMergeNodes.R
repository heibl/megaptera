## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2021-05-12)

#' @title Utilities for NCBI Taxdump
#' @description Insert a taxon as new node to the NCBI taxonomy.
#' @param x An object of class \code{\link{megapteraProj}} or a data frame in
#'   parent-child format.
#' @param nodes A vector of mode \code{"character"} giving the accepted names of
#'   the nodes to be merged.
#' @param accept A character string giving the accpeted taxon name of the new node.
#' @details \code{taxdumpMergeNodes} can be used in two ways: (1) If \code{x} is an
#'   object of class \code{"megapteraProj"}, the changes are made in the
#'   underlying postgreSQL database and \code{TRUE} or \code{FALSE} is returned
#'   indicating success or failure. (\bold{not yet implemented!}) (2) If \code{x} is a data frame in
#'   parent-child format, the changes are made to this table, which is then
#'   returned.
#' @return See Details section.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpChildren}},
#'   \code{\link{taxdumpDropTip}}, \code{\link{taxdumpHigherRank}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
#' @export

taxdumpMergeNodes <- function(x, nodes, accept = nodes[1], new.child, new.child.rank){

  #########################################
  ## MODE 1: changes in postgreSQL database
  #########################################
  if (inherits(x, "megapteraProj")){
    
    stop("implement me!")
    # if (!missing(tab)){
    #   id <- which(tab$rank == rank)
    #   taxon <- tab$taxon[id]
    #   parent <- tab$taxon[(id + 1):(nrow(tab))]
    #   names(parent) <- tab$rank[(id + 1):(nrow(tab))]
    # }
    # 
    # conn <- dbconnect(x)
    # 
    # ## Check if taxon is really missing
    # ## --------------------------------
    # SQL <- paste("SELECT * FROM taxonomy",
    #              "WHERE", wrapSQL(taxon, "taxon", "="))
    # SQL <- dbGetQuery(conn, SQL)
    # if (nrow(SQL)) {
    #   dbDisconnect(conn)
    #   message("'", taxon, "' already present in taxonomy")
    #   return(TRUE)
    # }
    # 
    # ## Retrieve parent data. This can be done over a vector of
    # ## ranks that are ordered from lower to higher (i.e. towards the
    # ## root)
    # ## The parent of lowest available rank is then used as an anchor point
    # ## -------------------------------------------------------------------
    # SQL <- paste("SELECT * FROM taxonomy",
    #              "WHERE", wrapSQL(parent, "taxon", "=", NULL))
    # SQL <- lapply(SQL, dbGetQuery, conn = conn)
    # id <- sapply(SQL, nrow)
    # if (sum(id) == 0) {
    #   message("parent taxa not present in taxonomy")
    #   return(FALSE)
    # }
    # 
    # ## Construct new lineage to insert
    # ## -------------------------------
    # ID <- "SELECT max(parent_id), max(id) FROM taxonomy"
    # ID <- dbGetQuery(conn, ID)
    # ID <- max(ID)
    # anchor <- which(id == 1)[1]
    # if (anchor == 1){
    #   ## Case 1: taxon can be inserted directly, e.g. genus is already present
    #   ## ---------------------------------------------------------------------
    #   obj <- data.frame(id = ID + 1, parent_id = SQL[[anchor]]$id, 
    #                     rank, taxon, status = "scientific name", stringsAsFactors = FALSE)
    # } else {
    #   ## Case 2: taxon cannot be inserted directly, e.g. genus is absent
    #   ## ---------------------------------------------------------------
    #   ID <- rev((1:anchor) + ID)
    #   obj <- data.frame(id = ID, 
    #                     parent_id = c(ID[-1], SQL[[anchor]]$id), 
    #                     rank = tab$rank[1:anchor], 
    #                     taxon = tab$taxon[1:anchor], 
    #                     status = "scientific name",
    #                     stringsAsFactors = FALSE)
    # }
    # 
    # ## Insert new node into database table
    # ## -----------------------------------
    # SQL <- paste(wrapSQL(obj$id, term = NULL, boolean = NULL),
    #              wrapSQL(obj$parent_id, term = NULL, boolean = NULL),
    #              wrapSQL(obj$rank, term = NULL, boolean = NULL),
    #              wrapSQL(obj$taxon, term = NULL, boolean = NULL),
    #              wrapSQL(obj$status, term = NULL, boolean = NULL),
    #              sep = ", ")
    # SQL <- paste("INSERT INTO taxonomy", 
    #              "(id, parent_id, rank, taxon, status)",  
    #              "VALUES (", SQL, ")")
    # lapply(SQL, dbSendQuery, conn = conn)
    # 
    # 
    # dbDisconnect(conn)
    # return(TRUE)
    
  } else {
    
    ############################################
    ## MODE 2: changes in parent-child table 'x'
    ############################################
    
    ## Check if taxon names are present in x
    ## -------------------------------------
    if (any(!nodes %in% x$taxon)){
      stop("implement me!")
    }
    accept <- match.arg(accept, nodes)
    syn <- nodes[!nodes %in% accept]
    if (length(syn) > 2){
      stop("implement me")
    }
    
    acc_id <- x$id[x$taxon == accept]
    syn_id <- x$id[x$taxon == syn]
    origin <- paste(Sys.info()["user"], "taxdumpMergeNodes", Sys.Date(), sep = ":")
    
    
    ## Graft children of 'syn' onto 'accept'
    ## -------------------------------------
    if (missing(new.child)){
      x$parent_id[x$parent_id == syn_id] <- acc_id
    } else {
      if (!new.child %in% x$taxon){
        x <- taxdumpAddNode(x, rank = new.child.rank, 
                            taxon = new.child, parent = accept,
                            origin = origin)
      }
      nc_id <- x$id[x$taxon == new.child]
      x$parent_id[x$parent_id == syn_id] <- nc_id
    }
    
    
    ## Merge 'syn' with 'accept'
    ## -------------------------
    x$status[x$id == syn_id] <- "synonym"
    if ("origin" %in% names(x)){
      x$origin[x$id == syn_id] <- origin
    }
    x$id[x$id == syn_id] <- acc_id
    
    ## Make sure 'accept' is accepted name
    ## -----------------------------------
    if (x$status[x$taxon == accept] != "scientific name"){
      stop("implement me")
    }
    
    cat(silver("Node '" %+% magenta$bold(syn) %+% "' merged as synonym into node '" %+% 
                 magenta$bold(accept) %+% "' [" %+% cyan(paste0("id:", acc_id)) %+% "] \n", 
               sep = ""))
    
    return(x)
  }
}