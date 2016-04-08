## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2015-02-10)

whereToInsert <- function(phy, tax, tip, tip.rank = "spec", quiet = FALSE){
  
  ## create a subset of tax taht contains the tips of the tree
  ## plus the tip to be added
  tax <- rbind(tax[tax[, tip.rank] == tip, ],
               tax[tax[, tip.rank] %in% phy$tip.label, ])
  ## delete ranks where tip as no information
  tax <- tax[ , tax[1, ] != "-"]
  id <- apply(tax, 2, function(x) x[1] %in% x[-1])
  id <- max(which(id))
  rank <- names(tax)[id]
  if ( !quiet ) cat(rank, levels(tax[1, id])[tax[1, id]])
  id <- which(tax[, id] == tax[1, id])
  clade <- tax[, tip.rank][id[-1]]
  an <- noi(phy, clade)
  if ( !quiet ) cat(" (node ", an, ") ", sep = "")
  
  ## check monophyly
  oops <- setdiff(descendants(phy, an, labels = TRUE), clade)
  if ( length(oops) > 0 & length(clade) > 1 ){
    if ( length(clade) / length(oops) > .5 ){
      if ( !quiet ) cat("- WARNING: ", rank, "is paraphyletic:", paste(head(oops, 3), collapse = ", "),
          "[", length(oops), "]")
    } else {
      if ( !quiet ) cat("- WARNING: ", rank, "seems polyphyletic:", paste(head(oops, 3), collapse = ", "),
          "[", length(oops), "]")
    }
  } else {
    if ( !quiet ) cat("- ", rank, "is monophyletic")
  }
  an
}


