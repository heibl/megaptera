\name{supermatrix}
\alias{supermatrix}
\title{Retrieval and Concatentation of Alignments}
\description{
Concatenate alignmenent of individual loci into a supermatrix.
}
\usage{
supermatrix(megapteraProj, min.n.seq = 3, masked = TRUE, 
            partition, subset.locus, coverage.locus = 0.5, 
            exclude.locus, exclude.species,
            core.locus, core.species, squeeze.outgroup)
}
\arguments{
  \item{megapteraProj}{
an object of class \code{\link{megapteraProj}}.
  }
  \item{min.n.seq}{
numeric, the minimum number of sequences in any alignment that is required to be included in the supermatrix.
  }
  \item{masked}{
logical, if \code{TRUE} masked (if available) alignments will be concatenated; see \code{\link{stepH}}.  
  }
  \item{partition}{
a named list, each element of which contains the names of the loci to group together in one partition. The default is to assign each locus to its own partition. Partitioning of first/second and third nucleotides positions is not yet possible, but might be implemented in the future. 
  }
  \item{subset.locus}{
to be added.
  }
  \item{coverage.locus}{
numeric between 0 and 1 giving the required minimum coverage of any species in any alignment.
  }
  \item{exclude.locus}{
a vector of mode \code{"character"} giving the names of the loci in the database that will be excluded from the concatenation.
  }
  \item{exclude.species}{
\emph{currently unused.}
  }
  \item{core.locus}{
a vector of mode \code{"character"} giving the names of the 'core' loci: The resulting supermatrix will only contain species that are contained in the 'core' loci. This option is intended to create denser supermatrices.
  }
  \item{core.species}{
\emph{currently unused.}
  }
  \item{squeeze.outgroup}{
numeric, can be given to reduce the number of outgroup species: The function will select the \code{squeeze.outgroup} outgroup species with the best coverage. The idea is to create a more densely sampled outgroup.
  }
}
%\details{}
\value{

a list with four elements:
\item{zipname}{a character string giving the name of the files produced.}
\item{supermatrix}{a matrix of class \code{\link{DNAbin}}.}
\item{outgroup}{a vector of mode \code{"character"} giving the names of the species used as an outgroup.}
\item{partitions}{a vector of mode \code{"character"} giving the partitions of the supermatrix in the same spelling as accepted by RAxML.}

In addition three files are written to the working directory: (1) a NEXUS-formatted file and (2) a PHYLIP-formatted file containing the sequence matrix and (3) a zipped directory containing the PHYLIP-formatted sequence matrix plus output and partitions as separate ASCII-formatted files.
}
%\references{}
\author{Christoph Heibl}
%\seealso{}
%\examples{}
