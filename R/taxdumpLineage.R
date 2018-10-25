## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2018-06-13)

#' @title Utilities for NCBI Taxdump
#' @description Get all higher ranks including a given taxon.
#' @param tax Either an object of class \code{\link{megapteraProj}} or a data
#'   frame as returned by \code{\link{dbReadTaxonomy}}.
#' @param taxon A character string giving the name of the taxon.
#' @seealso \code{\link{taxdumpAddNode}}, \code{\link{taxdumpChildren}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}}
#' @export

taxdumpLineage <- function(tax, taxon){
  
  ## Get taxonomy if necessary (i.e. tax is not a parent-child-table)
  ## --------------------------------------------------------------
  if (inherits(tax, "megapteraProj")){
    tax <- dbReadTaxonomy(tax)
  }
  
  ## Determine which separator is used by 'tax' and impose it on 'taxon'
  ## -------------------------------------------------------------------
  ## Beware of evil strings like 'Tuberculina sp. Ru_hy-01'
  test <- head(tax$taxon[tax$rank == "species"])
  underscore <- length(grep("[[:upper:]][[:lower:]+]_[[:lower:]]", test)) > 0
  if (underscore){
    taxon <- gsub(" ", "_", taxon)
  } else {
    taxon <- gsub("_", " ", taxon)
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
  
  # while (!any(c("root", "Root of life") %in% obj$taxon)){
  while (!any(root %in% obj$taxon)){
    # pid <- tax[tax$id == pid$parent_id, c("parent_id", "id", "taxon", "rank")]
    pid <- tax[tax$id == pid$parent_id & tax$status == "scientific name", ]
    if (!nrow(pid)) stop(obj$taxon[nrow(obj)], " (id=", obj$id[nrow(obj)], ", pid=", 
                         obj$parent_id[nrow(obj)],") has no parent")
    obj <- rbind(obj, pid)
  }
  obj
}
