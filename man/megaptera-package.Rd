\name{megaptera-package}
\alias{megaptera-package}
\alias{megaptera}
\docType{package}
\title{
MEGA-Phylogenie TEchniques in R
}
\description{
This package automates the sequence retrieval and procession from internet databases (currently \bold{Nucleotide} at \bold{GenBank}) to produce (possibly large) aligned nucleotide sequence datasets.
}
\details{
\tabular{ll}{
Package: \tab megaptera\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
 Date: \tab 2018-01-11\cr
License: \tab GPL (>= 2)\cr
}
\bold{Welcome to megaptera!}
To work with \code{megaptera} you need a functional installation of \bold{postgreSQL} on your computer. Please see \url{http://www.postgresql.org} for details on installation and working with postgreSQL. Although not necessary to build datasets with \code{megaptera}, you might benefit from a basic knowledge of standard query language (SQL).
You also will need a program for \bold{multiple sequence alignment} (MSA). Currently, the only option here is MAFFT and I recommend you take a look at \url{http://align.bmr.kyushu-u.ac.jp/mafft/software/}.
In addition, poorly aligned parts of the dataset may be \bold{masked}. Although optional, this step is highly recommended. Again, there is currently one single program used by \code{megaptera}: Gblocks (\url{http://molevol.cmima.csic.es/castresana/Gblocks.html}).
If you are new to \code{megaptera}, I'd suggest to take a look the package vignette.
}
 
\author{
Christoph Heibl
Maintainer: Christoph Heibl <christoph.heibl@gmx.net>
}
\references{
Nucleotide at NCBI: \url{http://www.ncbi.nlm.nih.gov/nuccore}
}
\keyword{ package }
%\seealso{\code{\link[ape:ape-package]{ape}}}
%\examples{}
