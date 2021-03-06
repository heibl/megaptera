% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkSpecLocus.R
\name{checkSpecLocus}
\alias{checkSpecLocus}
\title{Barplot of Species Numbers per Locus}
\usage{
checkSpecLocus(
  megProj,
  stage = "selected",
  subset.taxon,
  subset.locus,
  plot = TRUE
)
}
\arguments{
\item{megProj}{An object of class \code{\link{megapteraProj}}.}

\item{stage}{A vector of mode \code{character}; valid are \code{"retrieved"},
\code{"sel"}, or \code{"blocks"} referring to the number of species (1)
found on in the database, (2) selected according to coverage and identity (see
\code{\link{locus}}); or included in the final alignment, respectively.}

\item{subset.taxon}{A vector of mode \code{"character"}, that can be used to choose
a subset of the total taxa available.}

\item{subset.locus}{A vector of mode \code{"character"}, that can be used to choose
a subset of the total loci available.}

\item{plot}{Logical, indicating if the barplot should be produced on the
current graphical device or not.}
}
\value{
A list, but mainly called for its side effect, the plotting of a
  barplot.
}
\description{
Produce a barplot showing the number of species found for each
  locus. Thereby the number of species either refers to the number of species
  found on GenBank, the number of species selected according to identity and
  coverage, or the number of species included in the final alignment.
  Additionally, black lines to right of the bars indicate the cumulative
  species numbers for each locus if this locus was to be concatenated with
  all other loci having more species.
}
\seealso{
\code{\link{checkMissingSpec}} to create a table of species that
  were missed/lost during consecutive steps of the pipeline.
}
