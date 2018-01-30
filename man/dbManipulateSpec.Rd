\name{dbManipulateSpec}
\alias{dbDeleteSpec}
\alias{dbExcludeSpec}
\alias{dbIncludeSpec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Manipulate Single Species in Database}
\description{
Manually manipulate problematic species that are not handled satisfactorily by the parameters \code{identity} and \code{coverage}: delete a species from the PostgreSQL database, exclude a species from alignment, or include (a formerly excluded) species in alignment.
}
\usage{
dbDeleteSpec(x, spec)

dbExcludeSpec(x, spec)

dbIncludeSpec(x, spec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{\link{megapteraProj}}.
}
  \item{spec}{
a character string giving a species name.
}
}
%\details{}
\value{None; as a side effect the table \code{"locus"} and \code{"spec_<gene>"} are updated.
}
% \references{}
\author{Christoph Heibl}
\note{
\code{dbExcludeSpec} must be called between \code{\link{stepF}} and \code{\link{stepMAFFT}}. Because \code{stepF} rewrites the \code{spec_<locus>} table, \code{dbExcludeSpec} must be called again after every use of \code{stepF}. Note that there is no warning if manual exclusions set with \code{dbExcludeSpec} are eliminated by \code{stepF}.
}
\seealso{\code{\link{dbUpdateTaxonomy}}, \code{\link{dbReadTaxonomy}},
\code{\link{dbWriteDNA}}, \code{\link{dbReadDNA}}
}
%\examples{}
