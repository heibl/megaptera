\name{cetacea}
\alias{cetacea}
\docType{data}
\title{Cetacea Phylogeny}
\description{
Phylogenetic tree of whales (Cetacea).
}
\usage{data("cetacea")}
\format{
  The format is:
List of 5
 $ edge       : num [1:178, 1:2] 91 92 93 94 95 96 97 98 99 100 ...
 $ Nnode      : int 89
 $ tip.label  : chr [1:90] "Delphinus_capensis" "Delphinus_delphis" "Tursiops_aduncus" "Stenella_frontalis" ...
 $ edge.length: num [1:178] 0.05338 0.00533 0.00864 0.0074 0.02409 ...
 $ node.label : chr [1:89] "" "100" "35" "57" ...
 - attr(*, "class")= chr "phylo"
 - attr(*, "order")= chr "cladewise"
}
%\details{}
\source{The phylogeny is based on a supermatrix of five loci assembled with \code{megaptera} on November 6, 2014 and was estimated under a GTR + Gamma model of sequence evolution in RAxML 8.0.20.}
%\references{}
\examples{
data(cetacea)
}
\keyword{datasets}
