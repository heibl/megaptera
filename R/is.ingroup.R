## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-03-28)

#' @title Does a Taxon Belong to Ingroup?
#' @description Check if taxon belong to the ingroup as defined by \code{\link{taxon}}.
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @param taxon A character vector containing taxon names.
#' @return A vector of logical values.
#' @export
#' @import DBI

is.ingroup <- function(megProj, taxon){
  
  ## CHECKS
  ## ------
  if (!inherits(megProj, "megapteraProj"))
    stop("'megProj' is not of class 'megapteraProj'")
  tip.rank <- megProj@taxon@tip.rank
  
  ingroup <- unlist(megProj@taxon@ingroup)
  if (!all(is.Linnean(ingroup))){
      ## ingroup genera compared against species names in stepD
      if (all(is.Linnean(taxon))){
        taxon <- strip.spec(taxon)
      }
      ingroup <- lapply(ingroup, taxdumpDaughters, x = megProj, tip.rank = tip.rank)
      ingroup <- lapply(ingroup, function(z, tr) z[z$rank == tr, "taxon"], 
                        tr = tip.rank)
      ingroup <- unlist(ingroup)
  }
  
  # outgroup <- unlist(megProj@taxon@outgroup)
  # if (!all(is.Linnean(outgroup))){
  #   stop("implement me!")
  # }
  
  taxon %in% ingroup
}
