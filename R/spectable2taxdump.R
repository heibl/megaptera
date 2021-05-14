## This code is part of the megaptera package 
##  © C. Heibl 2017 (last update 2021-05-11)

#' @title Utilities for NCBI Taxdump
#' @description Convert a list of species or genera into a taxonomy table in
#'   parent-child format.
#' @param x A vector of mode \code{"character"} giving a set of species or
#'   genus names.
#' @param root A character string giving the taxon name of the root node.
#' @param root.rank A character string giving the rank of the root node.
#' @param ranks A vector of mode \code{"character"} giving the ranks (= column
#'   names) that will be extracted from \code{x}. Ranks must be ordered from
#'   high to low (e.g. \code{c("family", "genus", "species")}).
#' @seealso \code{\link{dbReadTaxonomy}},
#'   \code{\link{taxdumpChildren}},\code{\link{taxdumpLineage}},
#'   \code{\link{taxdumpAddNode}}, \code{\link{taxdump2phylo}}.
#' @export

spectable2taxdump <- function(x, root = "root", root.rank = "no rank", ranks = "asis"){

  tax <- data.frame(parent_id = 1, 
                    id = 1:2, 
                    taxon = c("root", root), ## 'root' required by taxdumpLineage! 
                    rank = "no rank", 
                    status = "scientific name",
                    stringsAsFactors = FALSE)
  
  ## Create x from vector or list
  if (is.vector(x)){
    x <- data.frame(species = x, stringsAsFactors = FALSE)
  }
  
  ## Set ranks
  if (ranks == "asis"){
    ranks <- names(x)
  } else {
    ranks <- match.arg(ranks, names(x), several.ok = TRUE)
  }
  
  ## Genus can be derived from binomials
  if (!"genus" %in% ranks){
    x$genus <- strip.spec(x$species)
    ranks <- c(ranks, "genus")
    ranks[ranks %in% c("genus", "species")] <- ranks[match(c("genus", "species"), ranks)]
  }
  
  ## Handle no-binomial terms in species list (e.g. GenBank identifiers, etc)
  ## ------------------------------------------------------------------------
  id <- x$genus == x$species
  if (any(id)){
    x$species[id] <- paste("incertae-sedis", x$species[id], sep = "_")
    x$genus[id] <- "incertae-sedis"
  }

  ## Fill tax rank by rank
  ## ---------------------
  for (i in seq_along(ranks)){
    rr <- unique(x[[ranks[i]]])
    if (i == 1){
      pid <- 2
    } else {
      pid <- x[[ranks[i - 1]]][match(rr, x[[ranks[i]]])]
      pid <- tax$id[match(pid, tax$taxon)]
    }
    rr <- data.frame(parent_id = pid, 
                     id = seq_along(rr) + max(tax$id), 
                     taxon = rr, 
                     rank = ranks[i], 
                     status = "scientific name",
                     stringsAsFactors = FALSE)
    tax <- rbind(tax, rr)
  }
  
  ## Remove tag from no-binomial terms
  ## ---------------------------------
  tax$taxon[tax$rank == "species"] <- gsub("incertae-sedis_", "", tax$taxon[tax$rank == "species"])
  
  tax
}
