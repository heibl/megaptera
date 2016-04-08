\name{comprehensiveGuidetree}
\alias{comprehensiveGuidetree}
\title{Create Guide Tree for Alignment}
\description{
Create a complete guide tree for alignment of sequences based on the NCBI taxonomy togehter with a user-defined guide tree of higher ranks.
}
\usage{
comprehensiveGuidetree(megapteraProj, tip.rank, subset)
}
\arguments{
  \item{megapteraProj}{
an object of class \code{\link{megapteraProj}}.
  }
  \item{tip.rank}{
a character string giving the rank of the tips (e.g. \code{spec}, \code{gen}, ...)
  }
  \item{subset}{
a subset of species names the taxonomy should be limited to. Can be a DNA alignment of class \code{DNAbin}, a phylogenetic tree of class \code{phylo}, a data frame or a simple vector listing the species names.
  }
}
\details{It can be desirable to include a user-defined guide tree for two reasons. First, one might favor another phylogenetic hypothesis than inherent in the NCBI classification. Second, no classification can provide intra-rank phylogenetic relationships, like an explicit phylogeny does. Hence, via the the user-defined guide tree, it is possible to add much more resolution to the guide tree topology. 
}
\value{an object of class \code{\link{phylo}}.}
%\references{}
\author{Christoph Heibl}
%\note{}
\seealso{
\code{\link{dbReadTaxonomy}} for reading a taxonomic classification from the postgreSQL database.
}
%\examples{}
