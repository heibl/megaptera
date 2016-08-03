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
a subset of species names (Latin binomials) to which the taxonomy should be limited. Can be a DNA alignment of class \code{DNAbin}, a phylogenetic tree of class \code{phylo}, the name of a \code{spec_*} database table (as created by \code{\link{stepF}}), or a character vector.
  }
}
%\details{}
\value{For \code{dbReadTable} an object of class \code{"data.frame"} holding the taxonomic classification of target species.. \code{dbUpdateTaxonomy} is called only for its side effect creating and extending the \code{taxonomy} database table.}
%\references{}
\author{Christoph Heibl}
\seealso{\code{\link{ncbiTaxonomy}} for retrieval of taxonomies from the taxonomy database at NCBI; \code{\link{findRoot}} returns the lowest common rank for taxonomies.

\code{\link{dbReadDNA}} and \code{\link{dbWriteDNA}} to read and write DNA sequences. 


}

