## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-10-18)

#' @title Plot Very Large Phylogenies
#' @description Cut a phylogeny into suitable slices (at least try to - result
#'   will always depend on topoloy) and create a series of PDF files showing the
#'   individual slices. This technique is meant to facilitate the graphical
#'   inspection of very large phylogenies.
#' @param phy An object of class \code{\link{phylo}}.
#' @param file A vector of mode \code{"character"} giving a filename (and path)
#'   for the PDF file.
#' @return None, \code{slicePhylo} is called for its side effect of
#'   generating one or more PDF file.
#' @importFrom ips sister
#' @export

slicePhylo <- function(phy, file){
  
  i <- 1
  repeat {
    
    cat("\nFigure ", LETTERS[i], ": ", sep = "")
    
    ## compute largest monophyletic clade <= 750 sp
    ## --------------------------------------------
    n <- 1
    repeat{
      d <- sister(phy, n)
      dd <- union(n, d)
      if ( length(dd) > 750 ) break
      n <- dd
    }
    
    cat(length(n), "tips")
    
    ## {ips}noi does not work here because phy can contain
    ## duplicate tip labels, e.g. 'Genus p.p. (2/10 sp.)'
    nn <- n
    repeat {
      nn <- sort(unique(phy$edge[phy$edge[, 2] %in% nn, 1]))
      gg <- descendants(phy, min(nn))
      if ( all(n %in% gg) ) break
    } 
    pt <- extract.clade(phy, min(nn))
    pdfPhyloA0(pt, view = FALSE, save = TRUE, 
               file = paste0(file, "figure", LETTERS[i], ".pdf"))
    phy$tip.label[1] <- paste("FIGURE", LETTERS[i])
    phy <- drop.tip(phy, n[-1])
    
    i <- i + 1
    
    if ( Ntip(phy) <= 750 ){
      pdfPhyloA0(phy, view = FALSE, save = TRUE, 
                 file = paste0(file, "figure", LETTERS[i], ".pdf"))
      break
    }
  }
}