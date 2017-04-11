\name{stepH}
\alias{stepH}
\title{
Step H: Detect and Separate Unalignable Blocks 
}
\description{
Dependent on the substitution rate of the genomic region and the taxonomic depth of the study group, not always all accession can be aligned into a single alignment. \code{step H} will detect and separate such unalignable alignment blocks.
}
\usage{
stepH(x, max.mad)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  an object of class \code{\link{megapteraProj}}.
  }
  \item{max.mad}{
numeric, giving the treshold value for the assessment of saturation: alignments with a median average distance (MAD) of \code{max.mad} or greater will be broken into blocks. The default value of 0.01 has been estimated with simulation by Smith et al. (2009).
  }
}
%\details{}
\value{
None. \code{stepH} is called for its side effects%: (1) strings of DNA sequences with attribute data are stored in a pgSQL database; (2) a log file is written to the current working directory.
}
\references{
Smith, S. A., J. M. Beaulieu, and M. J. Donoghue. 2009. Mega-phylogeny approach for comparative biology: an alternative to supertree and supermatrix approaches. \emph{BMC Evolutionary Biology} \bold{9}:37.
}
\author{
Christoph Heibl
}
\seealso{
\code{\link{megapteraProj}}; \code{\link{stepG}} for the preceeding step and \code{\link{supermatrix}} for the concatenation of loci.
}
%\examples{}
