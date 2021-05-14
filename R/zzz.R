
#' @export

whichSplitsShared <- function(phy, cons){
  
  isSplitShared <- function(n, phy, cons){
    d <- descendants(phy, n, labels = TRUE)
    is.monophyletic(cons, d)
  }
  nn <- Ntip(phy) + (1:Nnode(phy))
  ## cannot test root node!
  c(TRUE, sapply(nn[-1], isSplitShared, phy = phy, cons = cons))
}
