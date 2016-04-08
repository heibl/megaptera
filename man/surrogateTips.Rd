\name{surrogateTips}
\alias{surrogateTips}
\title{
Congeneric Surrogate Species
}
\description{
Surrogate species in a phylogenetic tree are replaced by their congeneric study group members.
}
\usage{
surrogateTips(phy, species)
}
\arguments{
  \item{phy}{
 an object of class \code{\link{phylo}}.
}
  \item{species}{
 a character vector of species names. 
}
}
\details{
Sometimes not all of the species of interest are available for estimating a phylogenetic tree, either beacuse they have not been sequenced or because their sequences do not pass the sanity checks of \bold{megaptera} or another pipeline. In this case one can build the phylogeny using other species from the same genus than the missing species and use \code{surrogateTips} to replace the actual species names (the `surrogate') with the names of the missing species. This is done in a random way, such that repeated calls with create different phylogenetic trees and can be used to create a sample distribution of the surrogate clades' tree space.  

This are the assumptions that warrant the use of \code{surrogateTips}:
\enumerate{
\item The deepest split in the phylogeny must be considerably deeper than the diversification of genera.
\item The species that serve as surrogates must be monophyletic. (This is checked by the function.)
\item There should be little doubt that surrogate and focal species belong to the same genus.
}

}
\value{
an object of class \code{\link{phylo}}.
}
%\references{}
\author{Christoph Heibl}
\note{
Beware that the use of \code{surrogateTips} entails a rough tree manipulation. If you use this function for research make sure that you understand how it works and its implications on the conclusions that can be drawn out of such analyses. 
}

\seealso{
\code{\link{addTips}}
}
%\examples{}
