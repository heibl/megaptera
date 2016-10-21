\name{ncbiGenome}
\alias{ncbiGenome}
\title{
Find Organelle Genomes on NCBI GenBank
}
\description{
Finds the most representative organelle genomes of a given taxon for use as reference sequences.
}
\usage{
ncbiGenome(x, organelle, n = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{\link{megapteraProj}}.
  }
  \item{organelle}{
character, either \code{"mitochondrion"}, \code{"chloroplast"}, or any unambiguous abreviation of these. 
  }
  \item{n}{
numeric, the maximum number of genomes that will be chosen. Depending on the classification of the taxon as returned by \code{\link{stepA}} the   actual number of genomes returned can be less than \code{n}.
  }
}
\details{
\code{ncbiGenome} uses a four-step algorithm to produce a taxonomically balanced sample of reference organelle genomes: 

\enumerate{
\item Determine the root taxon for both ingroup and outgroup.

\item Find all organelle genomes present on NCBI GenBank for this taxa.

\item Using the taxonomic classifiaction, find the \code{n - x} basal lineages of the entire set of genomes; thereby \code{x} is often greater than 0 depending on the branching pattern (topology) encoded by the classifiaction.

\item For each lineage, randomly choose one organelle genome and return the results as data frame (see \code{Value} section).
}
}
\value{
a data frame with three columns:
\item{taxon}{scientific name as Latin binomial}
\item{gb}{UID: GenBank number}
\item{gi}{alternative UID (GIs will no longer be supported after august 2016!)}
}
\references{
NCBI Orgenelle Genome Resources:
\url{http://www.ncbi.nlm.nih.gov/genome/organelle/}
}
\author{Christoph Heibl}
\seealso{
\code{\link{locusRef}} to set reference sequences.
}
%\examples{}

