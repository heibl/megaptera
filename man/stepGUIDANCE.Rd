% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stepGUIDANCE.R
\name{stepGUIDANCE}
\alias{stepGUIDANCE}
\title{Detect Homology Uncertainty in Alingment}
\usage{
stepGUIDANCE(x, bootstrap = 100)
}
\arguments{
\item{x}{An object of class \code{\link{megapteraProj}}.}

\item{bootstrap}{An integer giving the number of alternative MSAs to be computed.}
}
\value{
None, \code{stepGUIDANCE} is called for its side effects.
}
\description{
Use GUIDANCE, GUIDANCE2 or HoT to calculate column-wise
  reliability scores for alignments.
}
\seealso{
\code{\link{megapteraProj}}; \code{\link{stepMAFFT}} for the
  preceeding step and \code{\link{stepH}} for the subsequent step.
}
