\name{stepBX}
\alias{stepBX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Step BX: Add External Sequences}
\description{
This step adds external sequences, i.e. sequences that were not retrieved from GenBank with \code{\link{stepB}}, to the postgreSQL database table \code{"acc_*"} (where * is the SQL-compatible name of the locus). Cannot be called before \code{stepB}.
}
\usage{
stepBX(x, dna, tag = "user-supplied", overwrite = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{\link{megapteraProj}}.
  }
  \item{dna}{
a list of class \code{\link{DNAbin}}.
  }
  \item{tag}{
a character string to tag the sequences in the database table.
  }
  \item{overwrite}{
logical, indicating if existing sequences in the database should be overwritten.
  }
}
\details{
It is important that the names of the sequences in \code{DNA} conform to the  pattern \code{"Genus_species_ID"}. The occurance of a underscore inside the strings \code{"Genus"}, \code{"species"}, and \code{"ID"} is illegal, but letters, numbers and hyphens are allowed.
}
\value{
None. \code{stepBX} is called for its side effects: (1) strings of DNA sequences with attribute data are stored in a pgSQL database; (2) a log file is written to the current working directory.
}
%\references{}
\author{
Christoph Heibl
}
\note{
The function is under development and the results should always be checked.
}
\seealso{
\code{\link{megapteraProj}}; \code{\link{stepB}} for the preceeding and \code{\link{stepC}} for the subsequent step.
}
%\examples{}
