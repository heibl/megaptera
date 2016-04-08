\name{megapteraPars-class}
\docType{class}
\alias{megapteraPars-class}
\title{Class "megapteraPars"}
\description{
Class for pipeline parameters of \code{\link{megapteraProj}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{megapteraPars(...)}.
}
%% Slots
\section{Slots}{
 \describe{
 \item{parallel}{
logical: if \code{TRUE}, several steps in the pipeline will be run in parallel, otherwise all steps are serial.
  }
  \item{cpus}{
numerical: if \code{TRUE}, several steps in the pipeline will be run in parallel, otherwise all steps are serial.
  }
  \item{cluster.type}{
character: if \code{TRUE}, several steps in the pipeline will be run in parallel, otherwise all steps are serial.
  }
  \item{update.seqs}{\emph{currently unused}.}
  \item{retmax}{numeric, giving the batch size when downloading sequences from the Entrez History server (default: 500).
}
  \item{max.gi.per.spec}{numeric, giving the maximum number of sequences that will be used per species. Can be used to avoid model organism (e.g., rice, \emph{Drosophila}, ...) cluttering up the pipeline with thousands of sequences (default: 100).
}
  \item{max.bp}{
numeric, the maximal length of DNA sequences in base pairs to be included in the alignment. The upper limit is determined by the alignment program and the specific alignment and can only be determined by trial-and-error (default: 5000).
}

  \item{reference.max.dist}{\emph{currently unused}.}

  \item{min.seqs.reference}{\emph{currently unused}.
  }
  \item{fract.miss}{
numeric, ranging between 0 and 1. To avoid long stretches of only a few sequences at the beginning and the ending of an alignment block a minimum required number of sequences can be set as a fraction of the total number of sequences in this alignment block. Has been superseeded by the \code{gb.*} parameters.
}

  \item{filter1}{\emph{currently unused}.}
  \item{filter2}{\emph{currently unused}.}
  \item{filter3}{\emph{currently unused}.}
  \item{filter4}{\emph{currently unused}.}

  \item{block.max.dist}{numeric, ranging between 0 and 1. \code{block.max.dist} gives the maximum genetic distance (measured as the fraction of divergent nucleotide positions) allowed in a sequence alignment block. The alignment of individual marker is iteratively broken into smaller blocks until this condition is met with.
  }
  \item{min.n.seq}{numeric, the minimum number of sequences required for an alignment block. Alignment blocks with less than \code{min.n.seq} are dropped from the output.
  }
  \item{gb1, gb2, gb3, gb4, gb5}{parameters for masking of alignment blocks with \code{\link{gblocks}}.
  }
%  \item{gb2}{parameter for cleaning of alignment blocks with \code{gblocks}}
%  \item{gb3}{}
%  \item{gb4}{}
%  \item{gb5}{}
 }
}
\author{Christoph Heibl}
\seealso{
\code{\link{megapteraProj}} for creating a \code{megaptera} project.
}
\examples{
megapteraPars()
}
