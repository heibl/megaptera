\name{strip.spec}
\alias{strip.spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Truncation of Species Names
}
\description{
Delete epithets form species names.
}
\usage{
strip.spec(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a character string
}
}
\details{
\code{strip.spec} tries to guess the separating character used; both \code{" "} and \code{"_"} are possible.
}
\value{
a character string
}
\author{
Christoph Heibl
}
%\seealso{\code{\link{help}}}
\examples{
strip.spec("Megaptera_novaengeliae")

## separating characters cannot be mixed!
spec <- c("Megaptera_novaengeliae", "Megaptera novaeangliae")
strip.spec(spec) # does not work
strip.spec(gsub("_", " ", spec)) # this works
strip.spec(gsub(" ", "_", spec)) # this works
}
