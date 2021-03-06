% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbReadMSA.R
\name{dbReadMSA}
\alias{dbReadMSA}
\title{Read Sequences from Database}
\usage{
dbReadMSA(
  x,
  locus,
  taxon,
  regex = TRUE,
  label = c("taxon", "acc"),
  confid.scores = "all",
  row.confid = 0,
  col.confid = 0,
  ignore.excluded = TRUE,
  blocks = "ignore"
)
}
\arguments{
\item{x}{An object of class \code{\link{megapteraProj}}.}

\item{locus}{A vector of mode \code{"character"} giving a locus name. This
argument is optional and, if specified, will override the locus definition
in \code{x}.}

\item{taxon}{A vector of mode \code{"character"} used to choose a subset of
vailable taxa. This can be either one or more taxon names or a regular
expression.}

\item{regex}{Logical: if \code{TRUE}, the string given via \code{taxon} will
be interpreted as a regular expression (see \code{\link{regex}}).}

\item{label}{A vector of mode \code{"character"} giving the attributes for
constructing the sequences labels, e.g. \code{"taxon"}, \code{"acc"}, or both.}

\item{confid.scores}{A vector of mode \code{"character"} specifying if and
how confidence scores will be returned. Use \code{"ignore"} to omit
confidence scores. Other options are \code{"all"}, \code{"row.means"}, and
\code{"col.means"}. Confidence scores are returned as an attribute and can
be accessed by \code{attr(obj, "cs")}.}

\item{row.confid}{A real number in the interval \code{[0, 1]} giving
the confidence threshold for alignment rows (i.e. taxa, sequences); only rows  scoring equal
or greater to \code{min.confid} will be selected.}

\item{col.confid}{A real number in the interval \code{[0, 1]} giving
the confidence threshold for alignment columns (i.e. nucletide positions); only rows
(i.e. taxa, sequences) or columns  scoring equal
or greater to \code{min.confid} will be selected.}

\item{ignore.excluded}{\emph{Currently unused}.}

\item{blocks}{\emph{Currently unused}.}
}
\value{
An object of class \code{\link{DNAbin}}.
}
\description{
Reads selected and assembled sequences (see \code{\link{stepF}}
  and subsequent steps) from the database, specifically from the table
  'sequence_selected'.
}
\seealso{
\code{\link{dbReadDNA}} and \code{\link{dbWriteDNA}}
}
