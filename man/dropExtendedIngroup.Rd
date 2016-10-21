\name{dropExtendedIngroup}
\alias{dropExtendedIngroup}
\title{
Extract Ingroup Phylogeny from `Extended Ingroup'
}
\description{
Extracts the ingroup phylogeny by dropping all tips that belong to the extended ingroup. Potentially, extended ingroup tips can be used as surrogate species. See Details section for explication of the concepts `extended ingroup' and `surrogate species'.
}
\usage{
dropExtendedIngroup(x, phy, surrogate = FALSE)
}
\arguments{
  \item{x}{
an object of class \code{\link[megapteraProj-class]{megapteraProj}}.
}
  \item{phy}{
an object of class \code{phylo}.
}
  \item{surrogate}{
logical, \emph{- to be described -}
}
}
\details{
\emph{Describe concepts of extended ingroup and surrogate species}.
}
\value{
an object of class \code{phylo}.
}
%\references{}
\author{Christoph Heibl}
%\seealso{\code{\link{help}}}
%\examples{}
