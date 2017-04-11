\name{dbDNA}
\alias{dbReadDNA}
\alias{dbWriteDNA}
\docType{data}
\title{Read and Write DNA Sequences}
\description{
Read and write DNA Sequences from/to PostgreSQL Database.
}
\usage{
dbReadDNA(x, tab.name, taxon, regex = TRUE, max.bp,
          min.identity, min.coverage, ignore.excluded = TRUE,
          subtree = FALSE, blocks = "split", masked = FALSE)

dbWriteDNA(conn, tab.name, dna, enforce.binomial = TRUE,
  status)
}
\arguments{
  \item{x}{
an object of class \code{\link{megapteraProj}}.
  }
  \item{conn}{
a connection object as produced by \code{\link[DBI]{dbConnect}}.
  }
  \item{tab.name}{
a character string giving the name of the database table (gene).
  }
  \item{taxon}{
a character string giving the name of the taxon.
  }
  \item{regex}{
logical: if \code{TRUE}, the string given via \code{taxon} will be interpreted as a regular expression (see \code{\link{regex}}).
  }
  \item{dna}{
an object of class \code{\link[ape]{DNAbin}} to write to the database.
  }
  \item{enforce.binomial}{
Logical
  }
  \item{max.bp}{
an integer, only sequences shorter equal or shorter than \code{max.bp} will be returned.
  }
  \item{min.identity}{
a real number between 0 and 1, only sequences with a fraction of at least \code{min.identity} nucleotides that are identical with the reference sequence will be returned.
  }
  \item{min.coverage}{
a real number between 0 and 1, only sequences with a fraction of at least \code{min.coverage} base pairs (compared to the reference sequence) will be returned.
  }
  \item{status}{
a character string to be written to the \emph{status} field of the PostgreSQL table.
  }
  \item{ignore.excluded}{
logical: if \code{FALSE}, \code{dbReadDNA} will also return sequences that are marked in the \emph{status} field as 'excluded', 'too long', or 'too distant'.
  }
  \item{blocks}{
a character string indicating how to handle alignment blocks: \code{"split"} causes blocks to be returned as elements of a list. \code{"concatenate"} means, blocks will be concatendated and returned as a single alignment. In order to get a (potentially unaligned) list of sequences, choose \code{"ignore"}.
  }
  \item{subtree}{
  logical.
  }
  \item{masked}{
logical: if alignments have been masked in \code{\link{stepH}}, setting \code{masked=TRUE} will render only the unmasked alignment positions.
  }
}

%\details{}
\author{Christoph Heibl}

