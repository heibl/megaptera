\name{addTips}
\alias{addTips}
\alias{addSingleTip}
\title{Add Missing Species to a Phylogenetic Tree}
\description{
Species missing in a phylogenetic tree will be added to that tree according to their taxonomic classification. 
}
\usage{
addSingleTip(phy, tip, tip.rank = "spec", tax, insert = "crown", 
             stem.edge = 0.5, quiet = FALSE)
             
addTips(phy, tax, tip.rank = "spec", tag = NULL, quiet = FALSE)
}
\arguments{
  \item{phy}{an object of class \code{\link{phylo}}.
  }
   \item{tip}{a character string giving the name of the species to be added.
}
\item{tip.rank}{a character string giving the ranks of the tips in \code{phy}; must be present as a column name in \code{tax}.
}
  \item{tax}{a data frame containing the taxonomy of the species in the tree plus the species that should be added to that tree.
  }
   \item{insert}{a character string indicating the positions where the species is to be inserted: \code{"crown"}, \code{"stem"}, \code{"randomly"}, or any unambigous abbreviation of these. This option will only have an effect if \code{phy} contains more than one congeneric of \code{tip}.
}
 \item{stem.edge}{a real number greater than 0 and smaller than 1, which gives the fraction of the terminal branch length that will be assigned to the branch subtending the newly created MRCA of \code{tip} and its single congeneric species in \code{phy}. Will have no effect if \code{phy} contains more than one congeneric of \code{tip}.
  }
  \item{quiet}{logical, indicating if screen output should be suppressed.}
  \item{tag}{\emph{to be added}}
}
%\details{}
\value{
an object of class \code{\link{phylo}}.
}
\author{Christoph Heibl}

\seealso{
\code{\link{dbReadTaxonomy}} for reading a taxonomy table from the postgreSQL database.
}
\examples{
data(cetacea)

tre <- addTips(cetacea$tree,
               cetacea$tax)
}
