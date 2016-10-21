## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-09-21)

harmonizePhyloData <- function(phy, data){
  
  not.in.data <- setdiff(phy$tip.label, data$species) 
  cat(length(not.in.data), "species not in \'data\'")
  not.in.tree <- setdiff(data$species, phy$tip.label) # not in data
  cat("\n", length(not.in.tree), " species not in \'phy\'", sep = "")
  
  phy <- drop.tip(phy, not.in.data)
  data <- data[!data$species %in% not.in.tree, ]
  data <- data[match(phy$tip.label, data$species), ]
  
  list(phylo = phy, data = data)
}



