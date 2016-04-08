\name{locusRef-class}
\docType{class}
\alias{locusRef-class}
\title{Class "locusRef"}
\description{
Class for genetic input parameters for \code{\link{megapteraProj}}.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{locusRef(..., not, search.fields = c("gene", "title"), reference.gi)}.
}
%% Slots
\section{Slots}{
 \describe{
 \item{kind}{
 
 }
 \item{sql}{
a vector of mode \code{character} holding a SQL-compatible transformation of the first alias for naming database tables.
}
  \item{aliases}{
a vector of mode \code{character} giving strings that should be seached for; the first element is taken to name the corresponding postgreSQL tables.
}
  \item{not}{
a vector of mode \code{character} giving strings that should be excluded from the search results; corresponds to the use of \code{NOT} in GenBank query.
}
  \item{search.fields}{
a vector of mode \code{character} setting the search fields (or attributes) of the Nucleotide database to be searched for the strings specified via the \dots argument.
  }
  \item{use.genomes}{
logical, if \code{TRUE} sequences of loci will be extracted from annotated genomes of chloroplasts and mitochondria. Due to the existance of different annotation styles, which might (still) be incompatible with \bold{megaptera}, this option can be used to turn the extraction of sequences from whole genomes off.
  }
  \item{align.method}{
a character string giving the alignment method in MAFFT. Available accuracy-oriented methods for less than 200 sequences are \code{"localpair"}, \code{"globalpair"}, and \code{"genafpair"}; \code{"retree 1"} and  \code{"retree 2"} are for speed-oriented alignment. The default is \code{"auto"}, which lets MAFFT choose an opproriate alignment method.
  }
  \item{min.identity}{
numeric between 0 and 1, giving the minimum proportion of nucleotides required to be identical with the reference sequence in order to be included in an alignment (default: 0.75)
  }
  \item{min.coverage}{
numeric between 0 and 1, giving the minimum proportion of nucleotide positions a sequence must have in common with the reference sequence in order to be included in an alignment (default: 0.5)
  }
  \item{reference}{
a list of class \code{DNAbin} that contains reference sequences from GenBank.
  }
 }
}
%% Extends
\section{Extends}{
Class \code{"locus"}, directly.
}
%%\details{}
%%\references{}
\author{Christoph Heibl}
%%\note{}
\seealso{
\code{\link{dbPars}}, \code{\link{taxon}}, and \code{\link{megapteraPars}} for defining of database parameters, taxa, and the pipeline's parameters, respectively, and \code{\link{megapteraProj}} for the bundling of input data.
}
%\examples{}
