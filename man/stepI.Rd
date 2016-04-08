\name{stepI}
\alias{stepI}
\title{
Step I: Detect Homology Uncertainty in Alingment
}
\description{
Use gblocks to detect positions of uncertain homology and remove them from the alignment.
}
\usage{
stepI(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{\link{megapteraProj}}.
  }
}
%\details{}
\value{
None. \code{stepI} is called for its side effects%: (1) strings of DNA sequences with attribute data are stored in a pgSQL database; (2) a log file is written to the current working directory.
}
%\references{}
\author{
Christoph Heibl
}
\seealso{
\code{\link{megapteraProj}}; \code{\link{stepH}} for the preceeding step and \code{\link{supermatrix}} for the concatenation of loci.
}
%\examples{}
