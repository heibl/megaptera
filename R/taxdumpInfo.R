## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2021-05-13)

#' @title Utilities for NCBI Taxdump
#' @description Summarize information for a taxon name
#' @param tax A data frame in parent-child format.
#' @param taxon A character string giving the name of the taxon to query.
#' @param lineage.depth Integer giving the maximum number of higher taxa when
#'   displaying the lineage leading back to the root.
#' @seealso \code{\link{taxdumpAddNode}}, \code{\link{taxdumpChildren}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}}
#' @importFrom crayon %+% bold magenta silver
#' @export

taxdumpInfo <- function(tax, taxon, lineage.depth){
  
  if (inherits(tax, "megapteraProj")){
    ## Need to set drop.extinct = FALSE, in order not to loose potential parentials
    ## that just have been added
    tax <- dbReadTaxonomy(tax, drop.extinct = FALSE)
    cat(silver("  [data read from database taxonomy table]\n"))
  } else {
    cat(silver("  [data read taxonomy table in workspace]\n"))
  }
  
  if (missing(taxon)){
    cat(silver("number of taxonomic concepts: " %+% magenta$bold(length(tax$id)) %+% "\n"))
    cat(silver("number of taxonom names:      ") %+% magenta$bold(length(tax$taxon)) %+% "\n")
    cat(silver("    of these accepted:        ") %+% 
          magenta$bold(nrow(tax[tax$status == "scientific name", ])) %+% "\n")
    cat(silver("    of these synonyms:        ") %+% 
          magenta$bold(nrow(tax[tax$status == "synonym", ])) %+% "\n")
    
    nof <- function(tax, rank){
      nrow(tax[tax$rank == rank & tax$status == "scientific name", ])
    }
    cat(silver("number of species:            ") %+% 
          magenta$bold(nof(tax, "species")) %+% "\n")
    cat(silver("number of genera:             ") %+% 
          magenta$bold(nof(tax, "genus")) %+% "\n")
    cat(silver("number of families:           ") %+% 
          magenta$bold(nof(tax, "family")) %+% "\n")
  
    return()
  }
  
  ## CHECKS
  if (length(taxon) > 1){
    taxon <- taxon[1]
    warning("'taxon' has more than one element - only the first element will be used")
  }
  hrule <- silver(paste0(paste(rep("-", 60), collapse = "") %+% "\n"))
  
  tt <- tax[tax$taxon %in% taxon, ]
  ss <- tax[tax$id == tt$id & tax$status == "synonym", ]
  ss <- paste(ss$taxon, collapse = ", ")
  ss <- ifelse(ss == "", "none", ss)
  
  cat(silver("name:         " ) %+% magenta$bold(tt$taxon) %+% "\n")
  cat(silver("rank:         ") %+% magenta$bold(tt$rank) %+% "\n")
  cat(silver("status:       ") %+% magenta$bold(tt$status) %+% "\n")
  cat(silver("origin:       ") %+% magenta$bold(tt$origin) %+% "\n")
  cat(silver("ID:           ") %+% magenta$bold(tt$id) %+% "\n")
  cat(silver("ID of parent: ") %+% magenta$bold(tt$parent_id) %+% "\n")
  cat(hrule)
  cat(silver("synonyms:     ") %+% magenta(ss) %+% "\n")
  cat(hrule)
  
  ## Children
  cc <- taxdumpChildren(tax, taxon, immediate = TRUE)
  cc <- sort(cc$taxon[!cc$taxon %in% taxon])
  if (!length(cc)) cc <- "none"
  for (i in seq_along(cc)){
    hh <- ifelse(i == 1, "children:     ", "              ")
    cat(silver(hh) %+% magenta(cc[i]) %+% "\n")
  }
  cat(hrule)
  
  ## Lineage
  cat(silver$bold("Lineage"))
  ll <- taxdumpLineage(tax, taxon)
  lll <- ll[-1, c("rank", "taxon")]
  if (!missing(lineage.depth)){
    lll <- head(lll, lineage.depth)
  }
  for (i in 1:nrow(lll)){
    indent <- paste0("\n", paste(rep(" ", i), collapse = ""))
    cat(indent %+% silver(paste0(lll$rank[i], ": ")) %+% magenta$bold(lll$taxon[i]))
  }
  cat("\n")
}
