nonmonophyletic <- function(phy, nodes){
  
  mrca <- noi(phy, nodes)
  ext <- descendants(phy, mrca)
  ext <- setdiff(ext, nodes)
  
  ext
}