## This code is part of the megaptera package
## © C. Heibl 2014 (2019-04-22)

#' @title Add Tips to a Phylogenetic Tree
#' @description Add tips (i.e. species) to a phylogenetic tree according to
#'   their taxonomic classification.
#' @param phy An object of class \code{\link{phylo}}.
#' @param tips A vector of mode \code{"character"} giving the names of the
#'   species to be added.
#' @param sister A vector of mode \code{"character"} giving the names of the
#'   species to which \code{tips} should be added as sister species. 
#' @param tax A data frame in parent-child format.
#' @param insert A character string indicating the positions where the species
#'   is to be inserted: \code{"crown"}, \code{"stem"}, \code{"randomly"}, or any
#'   unambigous abbreviation of these. This option will only have an effect for
#'   species that have in the phylogeny more than one sister taxon of the rank
#'   at which the are inserted.
#' @param ignore.monophyly Logical, indicating if monophyly should be considered
#'   when selecting the anchor point.
#' @param quiet Logical, indicating if screen output should be suppressed.
#' @details If only a species list is given with \code{tips}, the species are
#'   inserted based on the position of congeneric species in the phylogeney. As
#'   a consequence, species that have no congeners in the phylogeny cannot be
#'   inserted without giving information on the taxonomy of higher ranks, which
#'   can be done with \code{tax}.
#' @return An object of class \code{\link{phylo}}.
#' @seealso \code{\link{dbReadTaxonomy}} for reading a taxonomy table from the
#'   postgreSQL database.
#' @examples
#' data(cetacea)
#' missing_species <- setdiff(cetacea$tax$spec, cetacea$tree$tip.label)
#' tre <- addTips(phy = cetacea$tree, tips = missing_species, insert = "crown")
#' @export

addTips <- function(phy, tips, sister, tax, insert = "crown", 
                    ignore.monophyly = FALSE, quiet = FALSE){
  
  if (!inherits(phy, "phylo"))
    stop("'phy' is not of class 'phylo'")
  
  tips <- gsub(" ", "_", tips)
  
  if (missing(sister)){
    
    ## Case 1: Anchor point is taken from taxonomy
    ## -------------------------------------------
    if (missing(tax)){
      
      ## Species that have no congeners in the phylogeny can not be added without 
      ## taxonomic information of rank higher than 'genus'
      missing_genera <- setdiff(strip.spec(tips), strip.spec(phy$tip.label))
      id <- strip.spec(tips) %in% missing_genera
      not_addable <- tips[id]
      if (length(not_addable)){
        if (!quiet) cat("\n", length(not_addable), " species ", 
                        ifelse(length(not_addable) > 1, "have", "has"), 
                        " no congener in 'phy' and can not be added:",
                        paste("\n-", not_addable), sep = "")
        tips <- tips[!id]
      }
      
      ## Create taxononmy table from tips
      ## --------------------------------
      species <- union(tips, phy$tip.label)
      if (!quiet) cat("\nBuilding taxonomy table from species list")
      tax <- speclist2taxdump(species)
      
    } else {
      
      ## Check if 'tax' contains all 'tips'
      ## ----------------------------------
      tax$taxon <- gsub(" ", "_", tax$taxon)
      not <- setdiff(tips, tax$taxon)
      if (length(not)){
        if (!quiet) cat("\nWARNING: ", length(not), " elements of 'tips' not in 'tax'", 
                        paste0("\n- ", not), sep = "")
        tips <- setdiff(tips, not)
      }
      
      ## Create subset of taxonomy that fits phylogeny
      ## ----------------------------------------------
      tax <- taxdumpSubset(tax, species = union(phy$tip.label, tips))
    }
    if (!quiet) cat("\nNumber of species to add:", length(tips))
    
    ## Loop over species to add them to phylogeny; must be done like
    ## this because the phylogeny changes at each addition of species
    ## --------------------------------------------------------------
    for (i in tips){ 
      phy <- addSingleTip(phy, tip = i, insert = insert,
                          ignore.monophyly = ignore.monophyly,
                          tax = tax, quiet = quiet)
    }
  } else {
    ########################################################
    ## Case 2: Anchor point is sister species or sister genus
    ######################################################## 
    sister <- gsub(" ", "_", sister)
    if (length(tips) != length(sister)) stop("vectors 'tips' and 'sister' do not match")
    if (!quiet) cat("\nNumber of species to add:", length(tips))
    
    ## Loop over species to add them to phylogeny; must be done like
    ## this because the phylogeny changes at each addition of species
    ## --------------------------------------------------------------
    for (i in seq_along(tips)){ 
      # for (i in 1:12){ 
      # cat("\n", i)
      phy <- addSingleTip(phy, tip = tips[i], insert = insert,
                          ignore.monophyly = ignore.monophyly,
                          sister = sister[i], quiet = quiet)
    }
  }

  phy
}