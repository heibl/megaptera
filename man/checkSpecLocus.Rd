\name{checkSpecLocus}
\alias{checkSpecLocus}
\docType{data}
\title{
Barplot of Species Numbers per Locus
}
\description{
Produce a barplot showing the number of species found for each locus. Thereby the number of species either refers to the number of species found on GenBank, the number of species selected according to identity and coverage, or the number of species included in the final alignment. Additionally, black lines to right of the bars indicate the cumulative species numbers for each locus if this locus was to be concatenated with all other loci having more species.
}
\usage{checkSpecLocus(megapteraProj, stage = "sel", 
       plot = TRUE)}
\arguments{
  \item{megapteraProj}{
an object of class \code{\link{megapteraProj}}.
  }
  \item{stage}{
a vector of mode \code{character}; valid are \code{"gb"}, \code{"sel"}, or \code{"blocks"} referring to the number of species (1) found on GenBank, (2) selected according to coverage and identity (see \code{\link{locus}}); or included in the final alignment, respectively. 
  }
  \item{plot}{
logical, indicating if the barplot should be produced on the current graphical device or not.  
  }
}
%\details{}
\value{a list, but mainly called for its side effect, the plotting of a barplot.}
%\references{}
\author{Christoph Heibl}
\seealso{\code{\link{checkMissingSpec}} to create a table of species that were missed/lost during consecutive steps of the pipeline.}
%\examples{}

