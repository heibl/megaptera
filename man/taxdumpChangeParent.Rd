% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxdumpChangeParent.R
\name{taxdumpChangeParent}
\alias{taxdumpChangeParent}
\title{Utilities for NCBI Taxdump}
\usage{
taxdumpChangeParent(x, taxon, new.parent, orphaned.parent = "synonym")
}
\arguments{
\item{x}{An object of class \code{\link{megapteraProj}} or a
\code{data.frame} a parent-child-format.}

\item{taxon}{A character string giving the name of the taxon.}

\item{new.parent}{A character string giving the name of the new parent taxon.}

\item{orphaned.parent}{A character string indicating how to treat previous
parent taxa, that do not have other children than \code{"taxon"}, i.e.
which are orphaned.}
}
\value{
A \code{data.frame} in parent-child-format, if \code{class(x) ==
  "data.frame"}, or \code{TRUE}, if \code{class(x) == "megapteraProj"}. In
  the latter case \code{taxdumpChangeParent} is called for its side effect on
  the taxonomy table of the database.
}
\description{
Change a node's parent in the NCBI taxonomy.
}
\details{
Two choices are possible for \code{orphaned.parent}: (1)
  \code{"synonym"} will classify the previous, orphaned parent as a synonym
  of \code{"new.parent"}, whereas (2) \code{"delete"} will simply delete the
  orphaned parent from the taxonomy.
}
\seealso{
\code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy and
  \code{\link{dbReadTaxonomy}} for reading the project taxonomy; other
  taxdump-related tools: \code{\link{taxdumpChildren}},
  \code{\link{taxdumpDropTip}}, \code{\link{taxdumpHigherRank}},
  \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
  \code{\link{taxdump2phylo}} and \code{\link{taxdump_isTerminal}}.
}
