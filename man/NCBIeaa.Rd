\name{NCBIeaa}
\alias{NCBIeaa}
\docType{data}
\title{
NCBI Extended Amino Acid
}
\description{
This table provides a key to the Extended IUPAC Encoding of Amino Acids defined by NCBI.
}
\usage{data(NCBIeaa)}
\format{
  A data frame with 28 observations on the following 4 variables:
  \describe{
    \item{\code{Value}}{NCBI code}
    \item{\code{Symbol}}{One-letter abreviation of amino acids}
    \item{\code{Symbol3}}{Three-letter abreviation of amino acids}
    \item{\code{Name}}{Amino acids}
  }
}
\details{
The official IUPAC amino acid code has some limitations. One is the lack of symbols for termination, gap, or selenocysteine. Such extensions to the IUPAC codes are also commonly used by sequence analysis software. NCBI has created such a code which is simply the IUPACaa code above extended with the additional symbols.
}
\source{
\url{http://www.ncbi.nlm.nih.gov/books/NBK7198/#ch_datamod.datamodel.seqloc}
}
\references{
\url{http://www.ncbi.nlm.nih.gov/books/NBK7198/#ch_datamod.datamodel.seqloc}
}
\examples{
data(NCBIeaa)
}
\keyword{datasets}
