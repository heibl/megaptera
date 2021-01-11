## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2019-11-13)

#' @title Utilities for NCBI Taxdump
#' @description Get all children of certain rank for a given taxon.
#' @param tax Either an object of class \code{\link{megapteraProj}} or a data
#'   frame as returned by \code{\link{dbReadTaxonomy}}.
#' @param taxon A character string giving the name of the taxon.
#' @param immediate Logical, if \code{TRUE}, only the immediate children of
#'   \code{taxon} will be returned. The default (\code{FALSE}) return all
#'   children.
#' @param query.rank A character string giving the name of a rank that
#'   \code{"taxon"} must belong to. The default (\code{"any"}) does not imply
#'   any restriction on the rank of the taxon queried.
#' @param tip.rank A character string giving the name a rank. This rank will be
#'   treated as tip rank, i.e. all taxa of lower rank will be dicarded. 
#' @param status A character string defining the status of the taxon names to be
#'   returned, e.g. \code{"scientific name"} will return only currently
#'   accepted names, while \code{"all"} will return synonyms in addition.
#' @param indet A vector of character strings containing regular expressions
#'   (see Examples).
#' @param quiet Logical, use \code{quiet = TRUE} to suppress warning messages.
#' @return A data frame.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
#'   \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
#'   taxdump-related tools: \code{\link{taxdumpAddNode}},
#'   \code{\link{taxdumpDropTip}}, \code{\link{taxdumpHigherRank}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
#' @examples
#' # The set of default regular expressions used to identify nonvalid species binomials
#' indet.strings(collapse = TRUE)
#' @export

taxdumpChildren <- function(tax, taxon, immediate = FALSE, 
                            query.rank = "any", tip.rank = "species", 
                            status = "scientific name", 
                            indet, quiet = FALSE){
  
  ## checks
  ## ------
  if (length(taxon) > 1) stop("cannot handle more than one taxon")
  
  ## get taxonomy if necessary (i.e. tax is not a parent-child-table)
  ## --------------------------------------------------------------
  if (inherits(tax, "megapteraProj")){
    tax <- dbReadTaxonomy(tax)
  }
  
  if (status == "scientific name"){
    tax <- tax[tax$status == "scientific name", ]
  }
  
  ## Determine which separator is used by 'tax' and impose it on 'taxon'
  ## -------------------------------------------------------------------
  ## Beware of evil strings like 'Tuberculina sp. Ru_hy-01'
  test <- head(tax$taxon[tax$rank == "species"])
  underscore <- length(grep("^[[:upper:]][[:lower:]]{1,}_[[:lower:]]", test)) > 0
  if (underscore){
    taxon <- gsub(" ", "_", taxon)
  } else {
    taxon <- gsub("_", " ", taxon)
  }
  
  ## Check if taxon is contained in tax
  ## ----------------------------------
  if (!taxon %in% tax$taxon){
    stop("taxon '", taxon, "' is not present in 'tax'", sep = "")
  }
  
  ## If 'taxon' is of rank 'tip.rank', simply return it
  ## --------------------------------------------------
  if (tax$rank[tax$taxon == taxon] == tip.rank){
    return(tax[tax$taxon == taxon, , drop = FALSE])
  }
  
  ## Get id of 'taxon', i.e. the root node
  ## -------------------------------------
  if (length(grep("[[:alpha:]]", taxon)) == 1){
    if (query.rank == "any"){
      id <- tax[tax$taxon == taxon, "id"]
    } else {
      id <- tax[tax$taxon == taxon & tax$rank == query.rank, "id"]
    }
    if (!length(id)) {
      if (!quiet) warning("taxon '", taxon, "' not available")
      return(NULL)
    }
  } else {
    id <- taxon
  }
  
  ## Get all daughter nodes
  ## ----------------------
  this.id <- id
  gain <- length(id)
  while (gain > 0){
    this.id <- tax[tax$parent_id %in% this.id, "id"]
    id <- c(id, this.id)
    gain <- length(this.id)
    if (immediate) break
  }
  tax <- tax[tax$id %in% id, c(2, 1, 4, 3, 5)]
  rownames(tax) <- NULL
  
  ## Remove rows of rank below 'tip.rank'
  ## ------------------------------------
  id <- vector()
  this.id <- tax[tax$rank == tip.rank, "id"]
  gain <- length(this.id)
  while (gain > 0){
    this.id <- tax[tax$parent_id %in% this.id, "id"]
    id <- c(id, this.id)
    gain <- length(this.id)
  }
  tax <- tax[!tax$id %in% id, ]
  
  ## Remove *lineages* that do not terminate in structurally
  ## valid Latin binomials according to 'indet'
  ## ------------------------------------------
  if (missing(indet)) {
    indet <- indet.strings(collapse = TRUE)
  }
  notvalid <- grep(indet, tax$taxon)
  notvalid <- intersect(notvalid, which(tax$rank == "species"))
  message(notvalid)
  if (length(notvalid)){
    
    notvalid <- sort(tax$id[notvalid])
    for (i in notvalid){
      tax <- taxdumpDropTip(tax, i)
    }
    
    ## Might need this bit in the future if
    ## the above code turns out to be to slow
    ## --------------------------------------
    # pid <- tax$parent_id[notvalid]
    # pid <- table(pid)
    # z <- table(tax$parent_id)
    # z <- z[names(z) %in% names(pid)]
    # d <- z - pid
    # entire_genus <- d[d == 0]
    # pp_genus <- names(d)[d > 0]
    # 
    # notvalid <- intersect(notvalid, which(tax$parent_id %in% pp_genus))
    # tax <- tax[-notvalid]
    # 
    # notvalid <- grep(indet.strings(collapse = TRUE), tax$taxon)
    # notvalid <- intersect(notvalid, which(tax$rank == "species"))
  
  }
  tax
}