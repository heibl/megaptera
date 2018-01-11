\name{taxon}
\alias{taxon}
\alias{taxonGuidetree}

\title{Definition of Species/Higher Taxa}
\description{
Creates a S4 class defining a taxonomic parameters for a \code{\link{megapteraProj}}.
}
\usage{
taxon(ingroup, extend.ingroup = FALSE, outgroup,
  extend.outgroup = FALSE, kingdom, hybrids = FALSE,
  tip.rank = "species", reference.rank = "auto")
      
taxonGuidetree(ingroup, extend.ingroup = FALSE, outgroup,
  extend.outgroup = FALSE, kingdom, hybrids = FALSE,
  tip.rank = "species", reference.rank = "auto", guide.tree)
}
\arguments{
  \item{ingroup}{
a list of vectors of mode \code{"character"}, giving species names (Linnean binomials) or names of higher taxa that define the focal group. Linnean binomials and taxon names of higher rank must not be mixed in the same list. If the elements of the list are vectors of length greater than one, the first element of the vector is taken as the accepted name and the remaining are treated as synonyms. Synonyms are used in the search for sequences, but any results are stored in the database under the accepted name.
  }
  \item{extend.ingroup}{
logical, if \code{TRUE} and \code{ingroup} is a set of Linnean binomials, the search will be extended to include all congenerics of the species in the \code{ingroup}.  
  }
  \item{outgroup}{
a list of vectors of mode \code{"character"}, giving species names (Linnean binomials) or names of higher taxa that define the outgroup. Linnean binomials and taxon names of higher rank must not be mixed in the same list. If the elements of the list are vectors of length greater than one, the first element of the vector is taken as the accepted name and the remaining are treated as synonyms. Synonyms are used in the search for sequences, but any results are stored in the database under the accepted name.
  }
  \item{extend.outgroup}{
logical, if \code{TRUE} and \code{outgroup} is a set of Linnean binomials, the search will be extended to include all congenerics of the species in the \code{outgroup}.  
  }
  \item{kingdom}{
a vector of mode \code{"character"}, currently one of  \code{"Fungi"}, \code{"Metazoa"}, or \code{"Viridiplantae"}. The \bold{International Code of Nomenclature for algae, fungi, and plants} (ICBN) and the \bold{International Code of Zoological Nomenclature} (ICZN) do not exclude the possibility that genera and species of plants/fungi/algae and animals share the same name (e.g., \emph{Prunella} \emph{vulgaris}). Therefore it is necessary to include the kingdom into the search term when querying GenBank. 
}
  \item{hybrids}{
logical: if \code{TRUE}, hybrids (as recognized by the regular expression \code{"^x_|_x_"}) will be excluded from the pipeline.
  }
  \item{tip.rank}{
a vector of mode \code{"character"}, giving the rank of the tips. Currently, \code{"spec"} and \code{"gen"} are allowed for species and genera.
  }
  \item{reference.rank}{
a vector of mode \code{"character"}, giving the name of the rank to be used to create subsets of the sequences to derive the reference sequence(s). The default (\code{"auto"}) commits the selection of the reference rank to the pipeline and in most cases you should be fine using this option.
}
  \item{guide.tree}{an object of class \code{\link{phylo}} giving a user guide tree, which will be used for profile alignment.}
}
%\details{}
\value{
an object of class \code{taxon} or \code{taxonGuidetree}
}
%\references{}
\author{
Christoph Heibl
}
\seealso{
\code{\link{dbPars}}, \code{\link{locus}}, and \code{\link{megapteraPars}} for defining of database parameters, loci, and the pipeline's parameters, respectively, and \code{\link{megapteraProj}} for the bundling of input data.
}
%\examples{}
