\name{prune.phylo.rank}
\alias{prune.phylo.rank}
% TITLE
\title{Prune Phylogenies to Higher Ranks}
% DESCRIPTION
\description{
Prune tips of a certain taxonomic rank from an object of class \code{\link[ape]{phylo}} and obtain a new \code{phylo} object whose tips are of a higher rank.
}
\usage{
prune.phylo.rank(phy, tax, rank = "gen", add.nb = TRUE, 
   quiet = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{an object of class \code{\link[ape]{phylo}}.
}
  \item{tax}{a data frame containing taxonomic information 
    for the tip labels in \code{phy}.
}
  \item{rank}{a character string giving the name of a 
    column (= taxonomic rank) in \code{tax} to which 
    \code{phy} will be pruned.
}
  \item{add.nb}{logical
}
  \item{quiet}{logical}
}
%\details{}
% VALUE
\value{
an object of class \code{\link[ape]{phylo}}.
}
\author{
Christoph Heibl
}
\seealso{
\code{\link{addTips}} and \code{\link{addSingleTip}} to add terminal nodes to a phylogeny.
}
%\examples{}
