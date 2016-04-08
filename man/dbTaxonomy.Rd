\name{Taxonomy}
\alias{dbUpdateTaxonomy}
\alias{dbReadTaxonomy}
\title{Read and Write Taxonomy Table}
\description{
Create, extend and read the \code{taxonomy} table of a pgSQL database.
}
\usage{
dbUpdateTaxonomy(megapteraProj, taxonomy, tag)

dbReadTaxonomy(megapteraProj, subset, tag)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{megapteraProj}{
an object of class \code{\link{megapteraProj}}
  }
  \item{taxonomy}{
a data frame containing a taxonomic classification of the study species.
  }
  \item{tag}{
a character string used to tag the rows of \code{taxonomy} in the  postgreSQL table; this argument is mostly for internal use. 
  }
 \item{subset}{
a subset of species names the taxonomy should be limited to. Can be a DNA alignment of class \code{DNAbin}, a phylogenetic tree of class \code{phylo}, a data frame or a simple vector listing the species names.
  }
}
%\details{}
\value{For \code{dbReadTable} an object of class \code{"data.frame"} holding the taxonomic classification of target species.. \code{dbUpdateTaxonomy} is called only for its side effect creating and extending the \code{taxonomy} database table.}
%\references{}
\author{Christoph Heibl}
\seealso{\code{\link{ncbiTaxonomy}} for retrieval of taxonomies from the taxonomy database at NCBI. \code{\link{dbReadDNA}} and \code{\link{dbWriteDNA}} to read and write DNA sequences.
}

