\name{dbPars}
\alias{dbPars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set Database Parameters}
\description{
Sets the connection parameters for a PostgreSQL database and, as a side effect, creates the database, if it does not yet exist.
}
\usage{
dbPars(host = "localhost", port = 5432, dbname, 
    user = "postgres", password)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{host}{a vector of mode \code{"character"}
}
  \item{port}{numeric
}
  \item{dbname}{
a vector of mode \code{"character"}
}
  \item{user}{
a vector of mode \code{"character"}
}
  \item{password}{
a vector of mode \code{"character"}
}
}
\details{
\bold{megaptera} stores all data internally in a PostgreSQL database. \code{dbPars} represents the first step in setting up the \bold{megaptera} pipeline. See \code{\link[DBI]{dbConnect}} and \code{\link{PostgreSQL}} for further details about the connection parameters and procedure.
}
\value{
an object of class \code{dbPars}
}
\references{See the documentation at the PostgreSQL Web site \url{http://www.postgresql.org} for details.}
\author{Christoph Heibl}
%%\note{}
\seealso{
\code{\link{class.dbPars}} for a description of the class \code{dbPars}; \code{\link{taxon}}, \code{\link{locus}}, and \code{\link{megapteraPars}} for defining of taxa, loci, and the pipeline's parameters, respectively; and \code{\link{megapteraProj}} for the bundling of input data.
}
%\examples{}
