% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxdumpChildren.R
\name{taxdumpChildren}
\alias{taxdumpChildren}
\title{Utilities for NCBI Taxdump}
\usage{
taxdumpChildren(
  tax,
  taxon,
  immediate = FALSE,
  query.rank = "any",
  tip.rank = "species",
  status = "scientific name",
  indet,
  quiet = FALSE
)
}
\arguments{
\item{tax}{Either an object of class \code{\link{megapteraProj}} or a data
frame as returned by \code{\link{dbReadTaxonomy}}.}

\item{taxon}{A character string giving the name of the taxon.}

\item{immediate}{Logical, if \code{TRUE}, only the immediate children of
\code{taxon} will be returned. The default (\code{FALSE}) return all
children.}

\item{query.rank}{A character string giving the name of a rank that
\code{"taxon"} must belong to. The default (\code{"any"}) does not imply
any restriction on the rank of the taxon queried.}

\item{tip.rank}{A character string giving the name a rank. This rank will be
treated as tip rank, i.e. all taxa of lower rank will be dicarded.}

\item{status}{A character string defining the status of the taxon names to be
returned, e.g. \code{"scientific name"} will return only currently
accepted names, while \code{"all"} will return synonyms in addition.}

\item{indet}{A vector of character strings containing regular expressions
(see Examples).}

\item{quiet}{Logical, use \code{quiet = TRUE} to suppress warning messages.}
}
\value{
A data frame.
}
\description{
Get all children of certain rank for a given taxon.
}
\examples{
# The set of default regular expressions used to identify nonvalid species binomials
indet.strings(collapse = TRUE)
}
\seealso{
\code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
  \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
  taxdump-related tools: \code{\link{taxdumpAddNode}},
  \code{\link{taxdumpDropTip}}, \code{\link{taxdumpHigherRank}},
  \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
  \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
}
