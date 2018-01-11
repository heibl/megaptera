## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-10-18)

#' @title Plot Large Phylogenies
#' @description Create a PDF file of a large phylogeny.
#' @param phy An object of class \code{\link{phylo}}.
#' @param file A vector of mode \code{"character"} giving a filename (and path)
#'   for the PDF file.
#' @param view Logical, if \code{TRUE}, the PDF will be opened in the default
#'   PDF viewer, but nothing is saved.
#' @param save Logical, if \code{TRUE}, the PDF is saved to \code{file}.
#' @return None, \code{slicePhylo} is called for its side effect of generating a
#'   PDF file.
#' @importFrom ape nodelabels
#' @importFrom graphics plot
#' @export

pdfPhyloA0 <- function(phy, file = "bigtree.pdf",
                       view = FALSE, save = TRUE){
  
  ## graphical parameters optimized for A0
  cex = .4
  if ( Ntip(phy) > 750 ){
    height <- 46.81; width <- 33.11 # DIN A0
  }
  if ( Ntip(phy) <= 750 & Ntip(phy) > 380 ){
    height <- 33.11; width <- 23.41 # DIN A1
  }
  if ( Ntip(phy) <= 380 & Ntip(phy) > 190 ){
    height <- 23.41; width <- 16.56 # DIN A2
  }
  if ( Ntip(phy) <= 190 ){
    height <- 16.56; width <- 11.70 # DIN A3
  }
  
  ## tip colors
  tcol <- rep("black", Ntip(phy))
  tcol[grep("monotypic", phy$tip.label)] <- "blue"
  tcol[grep("incl[.]", phy$tip.label)] <- "grey35"
  tcol[grep("p[.]p[.]_-_[[:lower:]]", phy$tip.label)] <- "orange"
  tcol[grep("p[.]p[.]_-_[[:digit:]]", phy$tip.label)] <- "red"
  pdf(file, height = height, width = width)
  plot(phy, no.margin = TRUE,
       edge.width = .25,
       tip.color = tcol,
       cex = cex,
       type = "phylo", 
       use.edge.length = FALSE)
  nodelabels(phy$node.label, cex = cex, adj = c(1.1, -.3), frame = "n", col = "red")
  #tiplabels(cex = cex, adj = c(-.25, .5), frame = "n", col = "red")
  dev.off()
  if ( view ) system(paste("open", file))  
  if ( !save ) unlink(file)
}