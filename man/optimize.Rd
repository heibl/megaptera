\name{optimize}
\alias{optimizeCoverage}
\alias{optimizeIdentity}
\title{
Optimal Coverage and Identity Thresholds
}
\description{
For a given locus, get threshold values for minimum coverage and minimum identity that allow to include at least one sequence of every species in the database.
}
\usage{
optimizeCoverage(x, discard = 0, quiet = FALSE, plot = FALSE)

optimizeIdentity(x, discard = 0, quiet = FALSE, plot = FALSE)
}
\arguments{
  \item{x}{
an object of class \code{\link{megapteraProj}}.
  }
  \item{discard}{
integer, giving the number of species to `sacrifice'; e.g. if \code{discard = 2}, the threshold value is chosen such that all but two species will be selected for the final alignment.
  }
  \item{quiet}{
logical, indicating if screen output should be suppressed.
  }
  \item{plot}{
logical:
  }
}
%\details{}
\value{a list with three elements:
\item{current}{current threshold value}
\item{suggested}{suggest threshold value to include all species}
\item{ranking}{maximum identity/coverage values for each species in ascending order}
}
%\references{}
\author{Christoph Heibl}
\seealso{
\code{\link{locus}}, \code{\link{locusRef}}
}
%\examples{}
