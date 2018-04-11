## This code is part of the megaptera package
## © C. Heibl 2014 (2018-04-11)

#' @title Add Tips to a Phylogenetic Tree
#' @description Add tips (species) to a phylogenetic tree according to their
#'   taxonomic classification.
#' @param phy An object of class \code{\link{phylo}}.
#' @param tax A data frame in parent-child format.
#' @param tips A character string giving the names of the species to be added.
#' @param insert A character string indicating the positions where the species
#'   is to be inserted: \code{"crown"}, \code{"stem"}, \code{"randomly"}, or any
#'   unambigous abbreviation of these. This option will only have an effect if
#'   \code{phy} contains more than one congeneric of \code{tip}.
#' @param ignore.monophyly Logical, indicating if monophyly should be considered
#'   when selecting the anchor point.
#' @param quiet Logical, indicating if screen output should be suppressed.
#' @return An object of class \code{\link{phylo}}.
#' @seealso \code{\link{dbReadTaxonomy}} for reading a taxonomy table from the
#'   postgreSQL database.
#' @examples
#' data(cetacea)
#' # tre <- addTips(cetacea$tree, cetacea$tax) ## transform to parent-child-format!
#' @export

addTips <- function(phy, tax, tips, insert, ignore.monophyly = FALSE, quiet = FALSE){
  
  if (!inherits(phy, "phylo"))
    stop("'phy' is not of class 'phylo'")
  
  ## Check if 'tax' contains all 'tips'
  ## ----------------------------------
  not <- setdiff(tips, tax$taxon)
  if (length(not)){
    if (!quiet) cat("\nWARNING: ", length(not), " elements of 'tips' not in 'tax'", 
                    paste0("\n- ", not), sep = "")
    tips <- setdiff(tips, not)
  }
  
  ## Create sub set of taxonomy that fits phylogeny
  ## ----------------------------------------------
  tax <- taxdumpSubset(tax, species = union(phy$tip.label, tips))
  
  if (!quiet) cat("\nNumber of species to add:", length(tips))
  
  ## Loop over species to add them to phylogeny; must be done like
  ## this because the phylogeny changes at each addition of species
  ## --------------------------------------------------------------
  for (i in tips){ 
    phy <- addSingleTip(phy, tip = i, insert = insert,
                        ignore.monophyly = ignore.monophyly,
                        tax = tax, quiet = quiet)
  }
  phy
}