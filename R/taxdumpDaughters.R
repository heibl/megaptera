## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-05-28)

#' @title Utilities for NCBI Taxdump
#' @description Get all children of certain rank for a given taxon.
#' @param x Either an object of class \code{\link{megapteraProj}} or a data
#'   frame as returned by \code{\link{dbReadTaxonomy}}.
#' @param taxon A character string giving the name of the taxon.
#' @param tip.rank A character string giving the name a rank. This rank will be 
#'   treated as tip rank, i.e. all taxa of lower rank will be dicarded.
#' @param indet A vector of character strings containing regular expressions 
#'   (see Examples).
#' @return A data frame.
#' @seealso \code{\link{taxdumpLineage}}, \code{\link{taxdumpMRCA}}, \code{\link{taxdump2phylo}}.
#' @examples 
#' # The set of default regular expressions used to identify nonvalid species binomials
#' indet.strings()
#' @export

taxdumpDaughters <- function(x, taxon, tip.rank, indet = indet.strings()){
  
  ## checks
  ## ------
  if (length(taxon) > 1) stop("cannot handle more than one taxon")
  
  ## get taxonomy if necessary (i.e. x is not a parent-child-table)
  ## --------------------------------------------------------------
  if (inherits(x, "megapteraProj")){
    x <- dbReadTaxonomy(x)
  }
  
  ## get id of 'taxon', i.e. the root node
  ## -------------------------------------
  if (length(grep("[[:alpha:]]", taxon)) == 1){
    id <- x[x$taxon == taxon, "id"]
    if (!length(id)) {
      warning("taxon '", taxon, "' not available")
      return(NULL)
    }
  } else {
    id <- taxon
  }
  
  ## get all daughter nodes
  ## ----------------------
  this.id <- id
  gain <- length(id)
  while (gain > 0){
    this.id <- x[x$parent_id %in% this.id, "id"]
    id <- c(id, this.id)
    gain <- length(this.id)
  }
  x <- x[x$id %in% id, c(2, 1, 4, 3)]
  rownames(x) <- NULL
  
  ## Remove rows of rank below 'tip.rank'
  ## ------------------------------------
  id <- vector()
  this.id <- x[x$rank == tip.rank, "id"]
  gain <- length(this.id)
  while (gain > 0){
    this.id <- x[x$parent_id %in% this.id, "id"]
    id <- c(id, this.id)
    gain <- length(this.id)
  }
  x <- x[!x$id %in% id, ]
  
  ## Remove *lineages* that do not terminate in structurally
  ## valid Latin binomials according to 'indet'
  ## ------------------------------------------
  notvalid <- grep(indet.strings(collapse = TRUE), x$taxon)
  notvalid <- intersect(notvalid, which(x$rank == "species"))
  if (length(notvalid)){
    
    notvalid <- sort(x$id[notvalid])
    for (i in notvalid){
      x <- taxdumpDropTip(x, i)
    }
    
    ## Might need this bit in the future if
    ## the above code turns out to be to slow
    ## --------------------------------------
    # pid <- x$parent_id[notvalid]
    # pid <- table(pid)
    # z <- table(x$parent_id)
    # z <- z[names(z) %in% names(pid)]
    # d <- z - pid
    # entire_genus <- d[d == 0]
    # pp_genus <- names(d)[d > 0]
    # 
    # notvalid <- intersect(notvalid, which(x$parent_id %in% pp_genus))
    # x <- x[-notvalid]
    # 
    # notvalid <- grep(indet.strings(collapse = TRUE), x$taxon)
    # notvalid <- intersect(notvalid, which(x$rank == "species"))
  
  }
  x
}