\name{stepB}
\alias{stepB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Step B: Search and Download Sequences}
\description{
For any given project (see \code{\link{megapteraProj}}), this step searches the Nucleotide database on GenBank, downloads all sequences and stores them in a postgreSQL database table.
}
\usage{
stepB(x, update.seqs = "no")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{\link{megapteraProj}}.
  }
  \item{update.seqs}{
a character string determining the behaviour of \code{stepB} when it is run more than one time on the same locus/taxon combination: \code{"no"} means that only sequences that have been made available on GenBank since the last execution \code{stepB} will be downloaded; \code{"all"} means the database table with all sequences will be removed and all sequences will be downloaded de novo. Naturely, this option has no effect if \code{stepB} is called for the first time.
  }
}
\details{
All accessions are stored under their species name as appearing in the \emph{organism} field at GenBank, but information about infrageneric ranks is stripped of the taxon names before they are stored in the database.
}
\value{
None. \code{stepB} is called for its side effects: (1) strings of DNA sequences with attribute data are stored in a pgSQL database; (2) a log file is written to the current working directory.
}
\references{
\url{http://www.ncbi.nlm.nih.gov/books/NBK25501/}
}
\author{
Christoph Heibl
}
\note{
The function is under development and the results should always be checked.
}
\seealso{
\code{\link{megapteraProj}}; \code{\link{stepA}} for the preceeding and \code{\link{stepC}} for the subsequent step; \code{\link{stepBX}} for the addition of external sequences to the database.
}
%\examples{}
