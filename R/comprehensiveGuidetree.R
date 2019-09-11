## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-04-17)

#' @title Comprehensive Guide Tree
#' @description Creates a complete (or comprehensive) guide tree for the
#'   alignment of sequences based on the underlying taxonomic
#'   classification. If available (see \code{\link{taxonGuidetree}}), a
#'   user-defined guide tree of higher ranks will be incorporated to enhance the
#'   resolution of the comprehensive guide tree.
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @param tip.rank A character string giving the rank of the tips (e.g. 
#'   \code{species}, \code{genus}, ...).
#' @param subset A subset of species names the taxonomy should be limited to. 
#'   This can be a DNA alignment of class \code{DNAbin}, a phylogenetic tree of class
#'   \code{phylo}, a data frame or a character vector listing the species names.
#' @details It can be desirable to include a user-defined guide tree for two 
#'   reasons. First, one might favor another phylogenetic hypothesis than 
#'   inherent in the NCBI classification. Second, no classification can provide 
#'   intra-rank phylogenetic relationships, but an explicit phylogeny does. 
#'   Hence, via the user-defined guide tree, it is possible to add more 
#'   resolution to the guide tree topology.
#' @return An object of class \code{\link{phylo}}.
#' @seealso \code{\link{dbReadTaxonomy}} for reading a taxonomic classification 
#'   from the postgreSQL database.
#' @importFrom ape bind.tree drop.tip read.tree
#' @export

comprehensiveGuidetree <- function(megProj, tip.rank, subset){
  
  ## Default tip rank is the project's tip rank
  ## ------------------------------------------
  if (missing(tip.rank)) tip.rank <- megProj@taxon@tip.rank
  
  ## Read taxonomy table
  ## -------------------
  if (missing(subset)){
    tax <- dbReadTaxonomy(megProj, tip.rank = tip.rank, root = "mrca")
  } else {
    tax <- dbReadTaxonomy(megProj, tip.rank = tip.rank, subset = subset, root = "mrca")
  }
  tax <- tax[tax$status == "scientific name", ]
  # setdiff(ingroup_accepted, tax$taxon[tax$rank == "species"])
  # setdiff(tax$taxon[tax$rank == "species"], ingroup_accepted)
  tax$taxon <- gsub(" ", "_", tax$taxon)
  
  ## get outgroup taxa
  ## -----------------
  og <- gsub(" ", "_", unlist(megProj@taxon@outgroup))
  if (tip.rank == "genus" & all(is.Linnean(og))) og <- strip.spec(og)
  og <- intersect(og, tax$taxon) ## actual subset of outgroup
  
  ## collapse incertae sedis nodes
  ## -----------------------------
  incsed <- grep("incertae sedis", tax$taxon)
  if (length(incsed)){
    for (i in incsed){
      tax$parent_id[tax$parent_id == tax$id[i]] <- tax$parent_id[i]
    }
    tax <- tax[-incsed, ]
  }
  
  if (inherits(megProj@taxon, "taxonGuidetree")){
    
  ## A: graft tips on user-defined guidetree
  ## ---------------------------------------
    gt <- megProj@taxon@guide.tree
    gt$edge.length <- NULL
    
    ## Check, if any tips of the guide tree are not part of
    ## the project's taxonomy
    ## ----------------------
    drop_from_guidetree <- setdiff(gt$tip.label, tax$taxon)
    if (length(drop_from_guidetree)){
      ## If MRCA of subset is of lower rank than tips of 
      ## guidetree, it cannot be used
      if (!length(intersect(tax$taxon, gt$tip.label))) {
        warning("guidetree and taxon subset do not overlap")
        return(taxdump2phylo(tax, tip.rank))
      } else {
        gt <- drop.tip(gt, drop_from_guidetree)
      }
    }
    
    ## Create subtrees that will be plotted onto the guide
    ## tree's tips
    ## -----------
    subtrees <- lapply(gt$tip.label, taxdumpChildren, tax = tax, tip.rank = tip.rank)
    subtrees <- lapply(subtrees, taxdump2phylo, tip.rank = tip.rank)
    names(subtrees) <- gt$tip.label

    single.gen <- sapply(subtrees, is.character)
    multi.gen <- which(!single.gen)
    single.gen <- which(single.gen)
    
    ## Check, if subtrees contain all focal species/genera
    ## This check could be done earlier, eg. after step A
    ## ---------------------------------------------------
    test_present <- c(unlist(subtrees[single.gen]),
                      unlist(lapply(subtrees[multi.gen], function(z) z$tip.label)))
    test_required <- tax$taxon[tax$rank == tip.rank]
    test_missing <- setdiff(test_required, test_present)
    test_missing <- setdiff(test_missing, og) ## do not consider outgroup as missing
    if (length(test_missing)) stop("there is no anchorage in user-defined guide tree for\n- ", 
                                   paste(test_missing, collapse = "\n -"))
    
    ## graft subtrees onto guide tree
    ## ------------------------------
    if (length(single.gen)){
      gt$tip.label[single.gen] <- unlist(subtrees[single.gen])
    }
    if (length(multi.gen)){
      if (length(multi.gen) == 1 & !length(single.gen)){
        gt <- subtrees[[1]]
      } else {
        for (i in names(subtrees)[multi.gen]){
          gt <- bind.tree(gt, subtrees[[i]], where = which(gt$tip.label == i)) 
        } 
      }
    }
    
    ## Add outgroup if nessesary
    ## -------------------------
    cond1 <- any(tax[tax$rank == tip.rank, "taxon"] %in% og)
    cond2 <- !all(og %in% gt$tip.label)
    if (cond1 & cond2){
      tt <- read.tree(text = "(outgroup, ingroup);")
      gt <- bind.tree(tt, gt, which(tt$tip.label == "ingroup"))
      if (length(og) == 1){
        gt$tip.label[gt$tip.label == "outgroup"] <- og
      } else {
        
        og_phylo <- taxdumpSubset(megProj, species = og, root = "mrca")
        og_phylo <- taxdump2phylo(og_phylo, tip.rank)
        gt <- bind.tree(gt, og_phylo, which(tt$tip.label == "outgroup"))
      }
    }
    
  } else {
    
    ## B: no user-defined guidetree given, turn
    ## classification into guidetree
    ## -----------------------------
    # gt <- findRoot(megProj, "both")
    # gt <- taxdumpChildren(tax, tail(gt$taxon, 1), tip.rank)
    # gt <- taxdump2phylo(gt, tip.rank)
    gt <- taxdump2phylo(tax, tip.rank)
  }
  gt
}
