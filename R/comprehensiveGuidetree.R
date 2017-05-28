## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-05-28)

#' @title Comprehensive Guide Tree
#' @description Creates a complete (or comprehensive) guide tree for the
#'   alignment of sequences based on the underlying taxonomic
#'   classification. If available (see \code{\link{taxonGuidetree}}), a
#'   user-defined guide tree of higher ranks will be incorporated to enhance the
#'   resolution of the comprehensive guide tree.
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @param tip.rank A character string giving the rank of the tips (e.g. 
#'   \code{spec}, \code{gen}, ...)
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
#' @importFrom ape bind.tree drop.tip
#' @export

comprehensiveGuidetree <- function(megProj, tip.rank, subset){
  
  ## Default tip rank is the project's tip rank
  ## ------------------------------------------
  if (missing(tip.rank)) tip.rank <- megProj@taxon@tip.rank
  
  ## read taxonomy table
  ## -------------------
  if (missing(subset)){
    tax <- dbReadTaxonomy(megProj, tip.rank = tip.rank, root = "mrca")
  } else {
    tax <- dbReadTaxonomy(megProj, tip.rank = tip.rank, subset = subset, root = "mrca")
  }
  
  if (inherits(megProj@taxon, "taxonGuidetree")){
    
  ## A: graft tips on user-defined guidetree
  ## ---------------------------------------
    gt <- megProj@taxon@guide.tree
    gt$edge.length <- NULL
    
    ## check if any tips of the guide tree are not part of
    ## the project's taxonomy
    ## ----------------------
    drop_from_guidetree <- setdiff(gt$tip.label, tax$taxon)
    if (length(drop_from_guidetree)){
      gt <- drop.tip(gt, drop_from_guidetree)
    }
    
    ## create subtrees that will be plotted onto the guide
    ## tree's tips
    ## -----------
    subtrees <- lapply(gt$tip.label, taxdumpDaughters, x = tax, tip.rank = tip.rank)
    subtrees <- lapply(subtrees, taxdump2phylo, tip.rank = tip.rank)
    names(subtrees) <- gt$tip.label
    
    single.gen <- sapply(subtrees, is.character)
    multi.gen <- which(!single.gen)
    single.gen <- which(single.gen)
    
    ## add ingroup clade if nessesary
    ## ------------------------------
    # ingroup <- megProj@taxon@ingroup
    # ingroup <- sapply(ingroup, function(z) z[1]) # strip synonoms
    # ingroup <- gsub(" ", "_", ingroup)
    # ig.rank.id <- unique(which(ingroup == tax, arr.ind = TRUE)[, "col"])
    # ig <- unique(tax[tax[, ig.rank.id] %in% ingroup, rank.id])
    # ig <- ig[!ig %in% gt$tip.label]
    # if ( length(ig) > 0 ){
    #   ig.phylo <- ifelse(length(ig) == 1, 
    #                      ig,
    #                      paste("(", paste(ig, collapse = ","), ")", sep = ""))
    #   ig.phylo <- paste("(", ig.phylo, ");", sep = "")
    #   ig.phylo <- read.tree(text = ig.phylo)
    #   gt <- bind.tree(gt, ig.phylo, "root")
    # }
    
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
    
    ## add outgroup if nessesary
    ## -------------------------
    og <- unlist(megProj@taxon@outgroup)
    cond1 <- any(tax[tax$rank == tip.rank, "taxon"] %in% og)
    cond2 <- !all(og %in% gt$tip.label)
    if (cond1 & cond2){
      tt <- read.tree(text = "(outgroup, ingroup);")
      gt <- bind.tree(tt, gt, which(tt$tip.label == "ingroup"))
      if (length(og) == 1){
        gt$tip.label[gt$tip.label == "outgroup"] <- og
      } else {
        og_phylo <- taxdumpSubset(tax, og)
        og_phylo <- taxdump2phylo(og_phylo, tip.rank)
        tt <- bind.tree(tt, og_phylo, which(tt$tip.label == "outgroup"))
      }
    }
    
  } else {
    
  ## B: no user-defined guidetree given, turn
  ## classification into guidetree
  ## -----------------------------
    # gt <- findRoot(megProj, "both")
    # gt <- taxdumpDaughters(tax, tail(gt$taxon, 1), tip.rank)
    # gt <- taxdump2phylo(gt, tip.rank)
    gt <- taxdump2phylo(tax, tip.rank)
    }
  gt$tip.label <- gsub(" ", "_", gt$tip.label)
  gt
}