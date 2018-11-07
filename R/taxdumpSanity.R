## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2018-11-07)

#' @title Sanity Check for Parent-Child Taxonomic Tables
#' @description Does several sanity checks for taxonomic tables in parent-child format.
#' @param tax A data frame in parent-child format.
#' @details \code{taxdumpSanity} checks, if
#' \enumerate{
#' \item synonyms have the same parent as their corresponding accepted name
#' \item accepted species names the corresponding genus name assigned
#' }
#' @return Logical
#' @export

taxdumpSanity <- function(tax){

  is_sane <- TRUE
  cat("Number of taxon names:", nrow(tax))
  id_set <- sort(unique(tax$id))
  cat("\nNumber of taxon concepts:", length(id_set))
  
  ## Synonyms are defined by having the same id as accepted taxa
  with_syn <- length(unique(tax$id[tax$status == "synonym"]))
  cat("\nNumber of taxon concepts with synonyms:", with_syn)
  
  ## 1. Do synonyms have same parent as accepted taxa?
  ## -------------------------------------------------
  parents <- tapply(tax$parent_id, tax$id, function(z) length(unique(z)))
  parents <- names(parents)[parents > 1]
  if (length(parents)){
    cat("\nFATAL:", length(parents), "taxon concepts have more than one parent:", 
        formatSpecList(parents, n.element = 6))
    is_sane <- FALSE
  }
  
  ## 2. Find accepted species linked to the wrong genus. This may happen,
  ##    when recombination is done
  ## -----------------------------
  accepted_species <- tax$taxon[tax$rank == "species" & tax$status == "scientific name"]
  genus <- sapply(accepted_species, taxdumpHigherRank, x = tax, rank = "genus")
  id <- genus == strip.spec(accepted_species)
  if (!all(id)){
    n <- length(accepted_species[!id])
    cat("\nFATAL:", n, "accepted species", ifelse(n == 1, "name is", "names are"), 
        "linked to a non-corresponding genus name:", 
        formatSpecList(accepted_species[!id]))
    is_sane <- FALSE
  }
  
  cat("\n")
  is_sane
}
