## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-04-15)

surrogateTips <- function(phy, species){
  
  ## species present in tree:
  present <- species %in% phy$tip.label
  cat("..", length(which(present)), "of", length(present), "species present in input tree ..")
  cat("\n..", length(which(!present)), "missing species from input tree ..")
  
  ## genera and species present or absent in tree
  pSpec <- species[present]
  aSpec <- species[!present]
  pGen <- strip.spec(pSpec)
  aGen <- strip.spec(aSpec)
  
#   sGen <- aGen[aGen %in% pGen] # absent genera do have some species
#   aGen <- aGen[!aGen %in% pGen] # absent genera do not have any species
  
  phyGen <- strip.spec(phy$tip.label)
  commonGen <- intersect(aGen, phyGen)
   
  t1 <- table(aGen[aGen %in% commonGen])
  t2 <- table(phyGen[phyGen %in% commonGen])
  
  id <- t2 - t1

  ## prepare list of id's species in phy
  ## -----------------------------------
  gen2spec <- function(spec, gen){
    pattern <- paste("^", gen, "_", sep = "")
    spec[grep(pattern, spec)]
  }
  phySpec <- lapply(names(id), gen2spec, spec = phy$tip.label)
  names(phySpec) <- names(id)
  
  ## test monophyly of genera: surrogates will be placed only in
  ## monophyletic genera
  ## -------------------
  mp <- sapply(phySpec, is.monophyletic, phy = phy)
  id <- id[mp]

  insertSpec <- lapply(names(id), gen2spec, spec = species)
  names(insertSpec) <- names(id)
  
  ## drop excess tips form tree
  ## --------------------------
  tips2drop <- vector()
  for ( i in names(id) ){
    
    fixed <- intersect(pSpec, phySpec[[i]])
    free <- setdiff(phySpec[[i]], pSpec)
    insert <- setdiff(insertSpec[[i]], pSpec)
    
    n <- length(free) - length(insert)
    if ( n < 0 ){
      insert <- sample(insert, length(insert) - abs(n))
    }
    if ( n > 0 ){
      drop <- sample(free, n)
      free <- free[!free %in% drop]
      tips2drop <- c(tips2drop, drop)
    }
    if ( length(insert) != length(free) ) stop()
    phy$tip.label[phy$tip.label %in% free] <- insert
    
  }
  phy <- drop.tip(phy, tips2drop)
  
  ## species present in tree:
  present <- species %in% phy$tip.label
  cat("\n..", length(which(present)), "of", length(present), "species present in output tree ..")
  cat("\n..", length(which(!present)), "missing species from output tree ..")
  
  phy <- drop.tip(phy, setdiff(phy$tip.label, species))
  phy
}