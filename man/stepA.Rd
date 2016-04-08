\name{stepA}
\alias{stepA}
\title{
Step A: Download a NCBI Taxonomy
}
\description{Download a taxonomic classification of the focal taxa from the NCBI taxonomy database and store it in the postgreSQL database.}
\usage{
stepA(x)
}
\arguments{
  \item{x}{an object of class \code{\link{megapteraProj}}.
}
}
%\details{}
\value{None. \code{stepB} is called for its side effects: (1) a taxonomic classifiaction is stored in the pgSQL database; (2) a log file is written to the current working directory.}
%\references{}
\author{Christoph Heibl}
\seealso{\code{\link{megapteraProj}}}
%\examples{}
