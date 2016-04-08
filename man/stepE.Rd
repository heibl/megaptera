\name{stepE}
\alias{stepE}

\title{
Step E: Calculate Each Accession's Distance from Benchmark
}
\description{
In \code{stepE} the genetic distance of each sequence from the corresponding benchmark sequence will be calculated and written to the postgreSQL database.
}
\usage{
stepE(x)
}
\arguments{
  \item{x}{an object of class \code{\link{megapteraProj}}.}
}
%\details{}
%\value{}
%\references{}
\author{
Christoph Heibl
}
\note{
The function is under development and the results should always be checked.
}
\seealso{
\code{\link{megapteraProj}}; \code{\link{stepD}} for the preceeding and \code{\link{stepF}} for the subsequent step.
}
%\examples{}