## This code is part of the megaptera package
## © C. Heibl 2014 (2017-06-08)

#' @title Add Tips to a Phylogenetic Tree
#' @description Add tips (species) to a phylogenetic tree according 
#' to their taxonomic classification.
#' @param phy An object of class \code{\link{phylo}}.
#' @param tax A character string giving the name of the 
#' species to be added.
#' @param tip.rank A character string giving the ranks of the 
#' tips in \code{phy}; must be present as a column name in 
#' \code{tax}.
#' @param tag \emph{To be added.}
#' @param quiet Logical, indicating if screen output should be 
#' suppressed.
#' @return An object of class \code{\link{phylo}}.
#' @seealso \code{\link{dbReadTaxonomy}} for reading a taxonomy table 
#' from the postgreSQL database.
#' @examples 
#' data(cetacea)
#' tre <- addTips(cetacea$tree, cetacea$tax)
#' @export

addTips <- function(phy, tax, tip.rank = "spec", tag = NULL, quiet = FALSE){
  
  if ( !inherits(phy, "phylo") )
    stop("'phy' is not of class 'phylo'")
  
  ## there must be no synonym and tag columns
  tax$synonym <- NULL
  tax$tag <- NULL
  
  ## data frame containing species missing from phy
  ## ----------------------------------------------
  add <- tax[!tax[, tip.rank] %in% phy$tip.label, ]
  if ( !quiet ) cat("\nNumber of species to add:", nrow(add))
  
  ## add via genus
  add.tips <- union(add[, tip.rank], NULL) # union returns character!
  for ( i in add.tips ){ 
    if ( !quiet ) cat("\n.. ", i)
    phy <- addSingleTip(phy, tip = i, tip.rank = tip.rank, 
                        tax = tax, quiet = quiet)
  }
  
  ## append tag to randomly added species
  if ( !is.null(tag) ){
    id <- phy$tip.label %in% add.tips
    phy$tip.label[id] <- paste(phy$tip.label[id], tag, sep = "_")
  }
  phy
}