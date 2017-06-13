## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-06-08)

#' @title Utilities for NCBI Taxdump
#' @description Get all higher ranks including a given taxon.
#' @param x Either an object of class \code{\link{megapteraProj}} or a data
#'   frame as returned by \code{\link{dbReadTaxonomy}}.
#' @param taxon A character string giving the name of the taxon.
#' @export

taxdumpLineage <- function(x, taxon){
  
  ## get taxonomy if necessary (i.e. x is not a parent-child-table)
  ## --------------------------------------------------------------
  if (inherits(x, "megapteraProj")){
    x <- dbReadTaxonomy(x)
  }
  
  obj <- data.frame(stringsAsFactors = FALSE)
  pid <- x[x$taxon == taxon, c("parent_id", "id", "taxon", "rank")]
  obj <- rbind(obj, pid)
  
  if (!nrow(obj)) return(NULL)

  ## Subgenera can have the same name as genera
  ## e.g. Tabanus subg. Tabanus
  if (nrow(pid) == 2 & "subgenus" %in% pid$rank){
    pid <- pid[pid$rank != "subgenus", ]
    obj <- obj[obj$rank != "subgenus", ]
  }
  
  while (!"root" %in% obj$taxon){
    pid <- x[x$id == pid$parent_id, c("parent_id", "id", "taxon", "rank")]
    if (!nrow(pid)) stop(obj$taxon[nrow(obj)], " (id=", obj$id[nrow(obj)], ", pid=", 
                         obj$parent_id[nrow(obj)],") has no parent")
    obj <- rbind(obj, pid)
  }
  obj
}