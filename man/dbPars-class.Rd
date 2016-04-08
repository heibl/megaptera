\name{dbPars-class}
\docType{class}
\alias{dbPars-class}
\alias{class.dbPars}

\title{Class "dbPars"}
\description{
Class for database parameters for \code{\link{megapteraProj}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{dbPars(host = "localhost", port = 5432, dbname, user = "postgres", password)}.
}
\section{Slots}{
 \describe{
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
}
\section{Methods}{
  \describe{
\item{show}{\code{signature(object = "taxon")}: prints database parameter setting}
}
}
%\details{}
%\value{an object of class \code{taxon}}
%\references{}
\author{Christoph Heibl}
\seealso{
\code{\link{dbPars}}, \code{\link{locus}}, and \code{\link{megapteraPars}} for defining of database parameters, loci, and the pipeline's parameters, respectively, and \code{\link{megapteraProj}} for the bundling of input data.
}
%\examples{}
\keyword{classes}
