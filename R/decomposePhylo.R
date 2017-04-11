## This code is part of the megaptera package
## © C. Heibl 2016 (2016-11-16)

#' @export

decomposePhylo <- function(phy, k = 200){
  
  if ( !is.binary.tree(phy) ) phy <- multi2di(phy, random = FALSE)
  
  rt <- Ntip(phy) + 1

  dsc <- phy$edge[phy$edge[, 1] == rt, 2]
  subtrees <- lapply(dsc, descendants, phy = phy)
  names(subtrees) <- dsc
  
  st.sizes <- sapply(subtrees, length)
  while ( any(st.sizes > k) ){

    id <- which(st.sizes > k)[1]

    dsc <- phy$edge[phy$edge[, 1] == names(id), 2]
    sbt <- lapply(dsc, descendants, phy = phy)
    names(sbt) <- dsc
    subtrees <- c(subtrees[-id], sbt)
    st.sizes <- sapply(subtrees, length)
  }
  subtrees <- lapply(subtrees, function(x, phy) phy$tip.label[x], phy = phy)
  st.size <- sapply(subtrees, length)
  st.label <- paste("st", 1:length(subtrees), sep = "-")
  names(subtrees) <- st.label
  z <- rep(st.label, st.size)
  z <- data.frame(subtree = z,
                  taxon = unlist(subtrees),
                  stringsAsFactors = FALSE)
  z <- z[match(phy$tip.label, z[, 2]), ]
  phy$tip.label <- z$subtree
  pruner <- function(phy, subtree){
    which(phy$tip.label %in% subtree)[-1]
  }
  w <- unlist(lapply(unique(z$subtree), pruner, phy = phy))
  p <- drop.tip(phy, w)
  list(guidetree = p,
       subtree.set = z)
}