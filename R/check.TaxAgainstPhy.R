## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-03-24)

#' @importFrom ape drop.tip is.monophyletic

check.TaxAgainstPhy <- function(phy, tax, rank = NULL){
  
  if ( !inherits(phy, "phylo") )
    stop("argument 'phy' is not of class 'phylo'")
  
  ## intersect tax and phy
  ## ---------------------
  tax <- tax[tax$spec %in% phy$tip.label, ]
  phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% tax$spec])
  
  x <- split(tax$spec, tax[rank])
  
  obj <- data.frame(taxon = names(x))
  obj <- cbind(obj, 
               topo = sapply(x, is.monophyletic, phy = phy))
  obj$topo[obj$topo] <- "mono"
  obj <- cbind(obj, size = sapply(x, length))
  obj
}