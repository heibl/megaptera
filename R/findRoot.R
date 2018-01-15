## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-11-06)

#' @title Lineage Down to the Root
#' @description Finds the lineage from one taxon, or the most recent common
#'   ancestor of several taxa, down to the root of the Tree of Life.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param what A character string, either \code{"ingroup"}, \code{"outgroup"}, or \code{"both"}.
#' @return A data frame with three columns:
#'  \item{id}{the unique identifier of the taxon}
#'  \item{taxon}{the scientific name of the taxon}
#'  \item{rank}{the rank of the taxon}
#' @importFrom methods slot
#' @export

findRoot <- function(x, what){
  
  what <- match.arg(what, c("ingroup", "outgroup", "both"))
  if (what == "both"){
    tax <- c(x@taxon@ingroup, x@taxon@outgroup)
  } else {
    tax <- unlist(slot(x@taxon, what))
  }
  
  if (all(is.Linnean(unlist(tax)))){
    
    ## read taxonomy from database
    tax <- dbReadTaxonomy(x, subset = tax)
    id <- all_ids <- 1
    id <- setdiff(tax[tax$parent_id == id, "id"], id)
    
    while (length(id) == 1){
      all_ids <- c(all_ids, id)
      id <- tax[tax$parent_id == id, "id"]
    }
    
    r <- tax[match(rev(all_ids), tax$id), c("parent_id", "id", "taxon", "rank")]
  } else {
    r <- dbReadTaxonomy(x)
    r <- lapply(tax, taxdumpLineage, tax = r)
    r <- r[!sapply(r, is.null)]
    
    
    if (length(r) > 1){
      
      ## find common set of nodes
      ## ------------------------
      rr <- lapply(r, function(x) x$id)
      for (i in 2:length(rr)){
        rr[[i]] <- intersect(rr[[i - 1]], rr[[i]])
      }
      
      ## find lowest rank in common set of nodes
      ## ---------------------------------------
      obj <- data.frame(common = rr[[i]], rank = NA)
      for (j in 1:nrow(obj)){
        obj[j, "rank"] <- mean(sapply(r, function(z, cc) which(z$id == cc), 
                                      cc = obj$common[j]))
      }
      r <- unique(do.call(rbind, r))
      mrca.id <- obj$common[order(obj$rank)]
      r <- r[match(mrca.id, r$id), c("id", "taxon", "rank")]
      
    } else {
      r <- r[[1]][, c("id", "taxon", "rank")]
    }
  }
  r
}