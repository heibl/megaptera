\name{strip.infraspec}
\alias{strip.infraspec}
\alias{strip.spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Truncation of Species Names
}
\description{
Delete infraspecific names and epithets form species names.
}
\usage{
strip.spec(x)

strip.infraspec(x, force.underscore = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a character string
}
  \item{force.underscore}{logical, if \code{TRUE}, white space will be replaced by underscore
}
}
\details{
\code{strip.spec} and \code{strip.infraspec} try to guess the separating character used; both \code{" "} and \code{"_"} are possible.
}
\value{
a character string
}
\author{
Christoph Heibl
}
%\seealso{\code{\link{help}}}
\examples{

(name <- "Vipera aspis aspis")
(name2 <- strip.infraspec(name))
(name3 <- strip.spec(name2))

names.set <- c("Vipera aspis aspis",
               "Vipera_aspis_atra")       
sapply(names.set, strip.infraspec)
sapply(names.set, strip.infraspec, force.underscore = TRUE)

}
