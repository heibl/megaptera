% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxdumpLineage.R
\name{taxdumpLineage}
\alias{taxdumpLineage}
\title{Utilities for NCBI Taxdump}
\usage{
taxdumpLineage(x, taxon)
}
\arguments{
\item{x}{Either an object of class \code{\link{megapteraProj}} or a data
frame as returned by \code{\link{dbReadTaxonomy}}.}

\item{taxon}{A character string giving the name of the taxon.}
}
\description{
Get all higher ranks including a given taxon.
}