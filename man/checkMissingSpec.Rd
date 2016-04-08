\name{checkMissingSpec}
\alias{checkMissingSpec}

\title{Check for Missing Species}
\description{
Create a table of species that failed to pass steps in the pipeline.
}
\usage{
checkMissingSpec(x, provenance = ".", tag = "ingroup")
}
\arguments{
  \item{x}{
an object of class \code{\link[dbPars-class]{dbPars}} or \code{\link[megapteraProj-class]{megapteraProj}}.
  }
  \item{provenance}{
a characters string as used in the column \code{provenance} the \code{acc_<locus>} tables in the postgreSQL database. Depending on your actual pipeline setup these may be \code{"ncbi"}, \code{"ncbi annotated"}, or \code{"bold"} as well as any user-defined strings.
  }
  \item{tag}{
a characters string as used in the column \code{tag} in the \code{taxonomy} table in the postgreSQL database. Depending on your choices for parameter values in \code{\link{taxon}} or \code{\link{taxonGuidetree}} these may be \code{"ingroup (NCBI)"}, \code{"extended ingroup (NCBI)"}, \code{"outgroup (NCBI)"}, or \code{"extended outgroup (NCBI)"} as well as any user-defined strings.
  }
}
%\details{}
\value{
a data frame with two columns:
\item{spec}{names of missing species}
\item{status}{why species are missing; might be (1) \code{not.on.genbank} (species for which no sequences were found on GenBank), (2)
\code{not.selected} (species that were not selected due to below-threshold identity or coverage), or (3) \code{not.aligned} (species that were lost during the alignment process; happens rarely)
}
}
%\references{}
\author{Christoph Heibl}
\seealso{\code{\link{checkSpecLocus}} to create a barplot of the species number per locus.}
%\examples{}

