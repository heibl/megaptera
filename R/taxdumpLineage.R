## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2019-03-05)

#' @title Utilities for NCBI Taxdump
#' @description Get all higher ranks including a given taxon.
#' @param tax Either an object of class \code{\link{megapteraProj}} or a data
#'   frame as returned by \code{\link{dbReadTaxonomy}}.
#' @param taxon A character string giving the name of the taxon.
#' @param highest.rank A character string giving the highest rank returned; if
#'   missing (default) the highest rank returned is the root.
#' @details If \code{taxon} is coded as a synonym, it will be replaced by the
#'   corresponding accepted name and a warning will be issued.
#' @return A data frame.
#' @seealso \code{\link{taxdumpAddNode}}, \code{\link{taxdumpChildren}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}}
#' @export

taxdumpLineage <- function(tax, taxon, highest.rank){
  
  ## Get taxonomy if necessary (i.e. tax is not a parent-child-table)
  ## --------------------------------------------------------------
  if (inherits(tax, "megapteraProj")){
    tax <- dbReadTaxonomy(tax)
  }
  
  ## If given, 'highest.rank' must be available
  ## ------------------------------------------
  if (!missing(highest.rank)) {
    if (!(any(tax$rank %in% highest.rank))) stop("'highest.rank' not available")
  }
  
  ## Determine which separator is used by 'tax' and impose it on 'taxon'
  ## -------------------------------------------------------------------
  ## Beware of evil strings like 'Tuberculina sp. Ru_hy-01'
  # test <- head(tax$taxon[tax$rank == "species"])
  # underscore <- length(grep("^[[:upper:]][[:lower:]]{1,}_[[:lower:]]", test)) > 0
  # if (underscore){
  #   taxon <- gsub(" ", "_", taxon)
  # } else {
  #   taxon <- gsub("_", " ", taxon)
  # }
  ## Try a variant, because taxon names 'Melyridae Dasytinae'
  ## --------------------------------------------------------
  if (gsub("_", " ", taxon) %in% tax$taxon){
    taxon <- gsub("_", " ", taxon)
  } else {
    if (!gsub(" ", "_", taxon) %in% tax$taxon){
      stop("taxon ='", taxon, "' not in 't", tax, "'", sep = "")
    } else {
      taxon <- gsub(" ", "_", taxon)
    }
  }
  
  ## Try to guess root (DIRTY!)
  ## --------------------------
  if (any(tax$taxon == "root")){
    root <- "root"
  } else {
    root <- setdiff(tax$parent_id, tax$id)
    root <- tax[tax$parent_id %in% root, ]
    root <- root$taxon[!root$rank %in% c("species", "genus")]
  }
  
  ## Prepare data frame to hold lineage
  ## ----------------------------------
  obj <- data.frame(stringsAsFactors = FALSE)
  # pid <- tax[tax$taxon == taxon, c("parent_id", "id", "taxon", "rank")]
  pid <- tax[tax$taxon == taxon, ]
  obj <- rbind(obj, pid)
  if (!nrow(obj)) return(NULL)

  ## Subgenera can have the same name as genera
  ## e.g. Tabanus subg. Tabanus
  ## --------------------------
  if (nrow(pid) == 2 & "subgenus" %in% pid$rank){
    pid <- pid[pid$rank != "subgenus", ]
    obj <- obj[obj$rank != "subgenus", ]
  }
  
  ## What, if taxon is a synonym?
  ## ----------------------------
  if (obj$status == "synonym"){
    oname <- obj$taxon
    nname <- tax$taxon[tax$id == obj$id & tax$status == "scientific name"]
    obj$taxon <- nname
    obj$status <- "scientific name"
    warning("synonym '", oname, "' replaced by accepted name ('", nname, "')")
  }
  
  i <- 1 ## control for unexitable loops
  while (!any(root %in% obj$taxon)){
    if (i > 100) stop("loop without exit")
    pid <- tax[tax$id == pid$parent_id & tax$status == "scientific name", ]
    if (!nrow(pid)) stop(obj$taxon[nrow(obj)], " (id=", obj$id[nrow(obj)], ", pid=", 
                         obj$parent_id[nrow(obj)],") has no parent")
    obj <- rbind(obj, pid)
    if (!missing(highest.rank)){
      if (highest.rank %in% obj$rank) break
    }
    i <- i + 1
  }
  obj
}
