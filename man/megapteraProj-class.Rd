\name{megapteraProj-class}
\docType{class}
\alias{megapteraProj-class}

\title{Class "megapteraProj"}
\description{
Class for holding all input data of the pipeline.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{megapteraProj(db, taxon, locus, align.exe, mask.exe, params = megapteraPars(), update = FALSE)}.
}
%% Slots
\section{Slots}{
\describe{
  \item{db}{an object of class \code{\link{dbPars}}.}
  \item{taxon}{an object of class \code{\link{taxon}}.}
  \item{locus}{an object of class \code{\link{locus}}.}
  \item{align.exe}{
a vector of mode \code{"character"}, giving name of the alignment program; currently only \bold{MAFFT} is allowed
}
  \item{mask.exe}{
a vector of mode \code{"character"}, giving name of the alignment masking program; currently only \bold{Gblocks} is allowed
}
  \item{params}{an object of class \code{\link{megapteraPars}}.}
  \item{update}{
logical: if \code{TRUE}, the pipeline's steps are executed as if called fro the first time, i.e., possibly overriding data and setting that have been previously achieved.
  }
 }
}
%\details{}
%\references{}
\author{Christoph Heibl}
\seealso{
\code{\link{dbPars}}, \code{\link{taxon}}, \code{\link{locus}}, and \code{\link{megapteraPars}} for defining of database parameters, taxa, loci, and the pipeline's parameters, respectively.
}
%\examples{}
