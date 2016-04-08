\name{checkSpecies}
\alias{checkSpecies}
\docType{data}
\title{Get Information about Species}
\description{
Returns information about the 'fate' of a single species along the pipeline.
}
\usage{checkSpecies(megapteraProj, spec)}
\arguments{
  \item{megapteraProj}{
an object of class \code{\link{megapteraProj}}.
  }
  \item{spec}{
a character string giving the name of a species.  
  }
}
%\details{}
\value{a data frame, as a side effect a summary message is printed on the screen.}
%\references{}
\author{Christoph Heibl}
\seealso{\code{\link{checkMissingSpec}} to create a list of species that were missed/lost during consecutive steps of the pipeline.}
%\examples{}

