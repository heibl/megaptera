## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-12-05)

#' @title Comprehensive Guide Tree
#' @description Creates a complete (or comprehensive) guide tree for the
#'   alignment of sequences based on the underlying taxonomic
#'   classification. If available (see \code{\link{taxonGuidetree}}), a
#'   user-defined guide tree of higher ranks will be incorporated to enhance the
#'   resolution of the comprehensive guide tree.
#' @param megapteraProj An object of class \code{\link{megapteraProj}}.
#' @param tip.rank A character string giving the rank of the tips (e.g. 
#'   \code{spec}, \code{gen}, ...)
#' @param subset A subset of species names the taxonomy should be limited to. 
#'   Can be a DNA alignment of class \code{DNAbin}, a phylogenetic tree of class
#'   \code{phylo}, a data frame or a simple vector listing the species names.
#' @details It can be desirable to include a user-defined guide tree for two 
#'   reasons. First, one might favor another phylogenetic hypothesis than 
#'   inherent in the NCBI classification. Second, no classification can provide 
#'   intra-rank phylogenetic relationships, but an explicit phylogeny does. 
#'   Hence, via the user-defined guide tree, it is possible to add more 
#'   resolution to the guide tree topology.
#' @return An object of class \code{\link{phylo}}.
#' @seealso \code{\link{dbReadTaxonomy}} for reading a taxonomic classification 
#'   from the postgreSQL database.
#' @export

comprehensiveGuidetree <- function(megapteraProj, tip.rank, subset){
  
  if ( missing(subset) ){
    tax <- dbReadTaxonomy(megapteraProj)
  } else {
    tax <- dbReadTaxonomy(megapteraProj, subset = subset)
  }
  
  ## adjust order: rightmost column must be the highest rank
  ## -------------------------------------------------------
  if ( all(c("spec", "gen") %in% names(tax)) ){
    if ( which(names(tax) %in% "spec") > which(names(tax) %in% "gen") )
      tax <- tax[, ncol(tax):1]
    if ( missing(tip.rank) ) tip.rank <- "spec"
  } else {
    if ( !all(c("fam", "gen") %in% names(tax)) ) stop("no columns spec, gen, fam")
    if ( which(names(tax) %in% "gen") > which(names(tax) %in% "fam") )
      tax <- tax[, ncol(tax):1]
    if ( missing(tip.rank) ) tip.rank <- "gen"
  }
  
  og <- unique(tax[grep("outgroup", tax$tag), ])
  
  ## drop columns 'synonym' and 'tag'
  og$synonym <- NULL
  og$tag <- NULL
  tax$synonym <- NULL
  tax$tag <- NULL
  
  if ( inherits(x@taxon, "taxonGuidetree") ){
    
  ## A: graft tips on user-defined guidetree
  ## ---------------------------------------
    gt <- megapteraProj@taxon@guide.tree
    gt$edge.length <- NULL
    
    ## get subtrees of rank <rank.id> = rank of tips in gt
    tip.rank.id <- which(names(tax) %in% tip.rank)
    rank.id <- sapply(names(tax), function(col, t1, t2) any(t1 %in% t2[, col]), 
                      t1 = gt$tip.label, t2 = tax)
    rank.id <- which(rank.id)
    if ( length(rank.id) > 1 ) stop("ambiguous ranks")
    subtrees <- split(tax[, (rank.id + 1):tip.rank.id ], f = tax[, rank.id])
    subtrees <- lapply(subtrees, unique)
    
    ## check if subtree roots are missing from guide tree ...
    not.gt <- setdiff(names(subtrees), gt$tip.label)
    og <- unique(og[, rank.id]) ## ... considering the outgroup
    not.gt <- setdiff(not.gt, og)
    if ( length(not.gt) > 0 ){
      stop(paste("subtree roots missing from gt:", paste(not.gt, collapse = ", ")))
    }
    
    single.gen <- which(sapply(subtrees, nrow) == 1)
    multi.gen <- which(sapply(subtrees, nrow) > 1)
    if ( length(multi.gen) > 0 ){
      subtrees[multi.gen] <- lapply(subtrees[multi.gen], tax2tree, tip.rank = tip.rank)
    }
    
    ## add ingroup clade if nessesary
    ## ------------------------------
    ingroup <- megapteraProj@taxon@ingroup
    ingroup <- sapply(ingroup, function(z) z[1]) # strip synonoms
    ingroup <- gsub(" ", "_", ingroup)
    ig.rank.id <- unique(which(ingroup == tax, arr.ind = TRUE)[, "col"])
    ig <- unique(tax[tax[, ig.rank.id] %in% ingroup, rank.id])
    ig <- ig[!ig %in% gt$tip.label]
    if ( length(ig) > 0 ){
      ig.phylo <- ifelse(length(ig) == 1, 
                         ig,
                         paste("(", paste(ig, collapse = ","), ")", sep = ""))
      ig.phylo <- paste("(", ig.phylo, ");", sep = "")
      ig.phylo <- read.tree(text = ig.phylo)
      gt <- bind.tree(gt, ig.phylo, "root")
    }
    
    ## add outgroup if nessesary
    ## -------------------------
    if ( !all(og %in% gt$tip.label) ){
      og.phylo <- ifelse(length(og) == 1, 
                         og,
                         paste("(", paste(og, collapse = ","), ")", sep = ""))
      og.phylo <- paste("(", og.phylo, ",ingroup);", sep = "")
      og.phylo <- read.tree(text = og.phylo)
      gt <- bind.tree(og.phylo, gt, which(og.phylo$tip.label == "ingroup"))
    }
    
    ## drop tips from guide tree that do not
    ## have a subtree to graft on
    ## --------------------------
    drop.id <- gt$tip.label[!gt$tip.label %in% names(subtrees)]
    if ( length(drop.id) > 0 ){
      gt <- drop.tip(gt, drop.id)
    }
    
    ## graft subtrees onto guide tree
    ## ------------------------------
    if ( length(single.gen) > 0 ){
      single.gen <- do.call(rbind, subtrees[single.gen])
      single.gen <- data.frame(from = rownames(single.gen), 
                               to = single.gen[, tip.rank],
                               stringsAsFactors = FALSE)
      gt$tip.label[match(single.gen$from, gt$tip.label)] <- single.gen$to
    }
    if ( length(multi.gen) > 0 ){
      for ( i in names(subtrees)[multi.gen] ){
        gt <- bind.tree(gt, subtrees[[i]], where = which(gt$tip.label == i)) 
      } 
    }
  } else {
    
  ## B: no user-defined guidetree given, turn
  ## classification into guidetree
  ## -----------------------------
    gt <- tax2tree(tax, tip.rank = tip.rank)
  }
  gt
}