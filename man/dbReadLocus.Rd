\name{dbReadLocus}
\alias{dbReadLocus}
\title{Taxon-Locus-Crosstable}
\description{
Create a data frame that contains for each taxon and locus the number of accessions found on GenBank/BOLD/etc. and if the taxon was selected for alignment.
}
\usage{
dbReadLocus(x, provenance = ".", tag)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{\link[megapteraProj-class]{megapteraProj}}.
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
a data frame with two types of colums
 \item{gb_<locus>}{gives the number of accessions per taxon and locus}
 \item{sel_<locus>}{states if a taxon was selected for final alignment}
}
\author{Christoph Heibl}

\seealso{
\code{\link{checkSpecLocus}}
}
%\examples{}

