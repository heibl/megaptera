% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/megapteraPars.R
\name{megapteraPars}
\alias{megapteraPars}
\title{Create an Object of Class "megapteraPars"}
\usage{
megapteraPars(...)
}
\arguments{
\item{...}{Arguments in \code{tag = value} form. The tags must come from the
names of the parameters described in the `Pipeline Parameters' section.}
}
\description{
S4 Class for parameters of a megaptera project pipeline, as
  stored in \code{\link{megapteraProj}}.
}
\details{
The pipeline's verbosity can be fine-tuned with \code{debug.level}:
  \tabular{ll}{ 0 \tab No progess and diagnostic messages\cr 1 \tab Messages
  on screen\cr 2 \tab Messages logged to file\cr 3 \tab Messages on screen
  and logged to file\cr 4 \tab Same as 3, in addition current data is saved
  as .rda object in case of a foreseeable error\cr 5 \tab Same as 4, in
  addition current data is always saved\cr }
}
\section{Pipeline Parameters}{
 
\describe{ 
\item{\code{data.path}}{A character string giving the path to the directory there all
  data and results will be stored (see \code{\link{megapteraInit}}).}
\item{\code{gb.seq.download}}{A character string defining how sequences should be
  downloaded from GenBank Nucleotide; can be \code{"eutils"} or \code{"ftp"}.}
\item{\code{debug.level}}{Numeric, a
  number between 0 and 5, determining the pipeline's verbosity (see Details). (default: 1)}
\item{\code{parallel}}{Logical: if
  \code{TRUE}, several steps in the pipeline will be run in parallel,
  otherwise all steps are serial.} \item{\code{cpus}}{Numerical: if
  \code{TRUE}, several steps in the pipeline will be run in parallel,
  otherwise all steps are serial.} \item{\code{cluster.type}}{A character
  string: if \code{TRUE}, several steps in the pipeline will be run in
  parallel, otherwise all steps are serial.}
  \item{\code{update.seqs}}{\emph{Currently unused}.}
  \item{\code{retmax}}{Numeric, giving the batch size when downloading
  sequences from the Entrez History server (default: 500).}
  \item{\code{max.gi.per.spec}}{Numeric, giving the maximum number of
  sequences that will be used per species. Can be used to avoid model
  organism (e.g., rice, \emph{Drosophila}, ...) cluttering up the pipeline
  with thousands of sequences (default: 1000).} \item{\code{max.bp}}{Numeric,
  the maximal length of DNA sequences in base pairs to be included in the
  alignment. The upper limit is determined by the alignment program and the
  specific alignment and can only be determined by trial-and-error (default:
  5000).} \item{\code{reference.max.dist}}{\emph{Currently unused}.}
  \item{\code{min.seqs.reference}}{\emph{Currently unused}.}
  \item{\code{fract.miss}}{Numeric, ranging between 0 and 1. To avoid long
  stretches of only a few sequences at the beginning and the ending of an
  alignment block a minimum required number of sequences can be set as a
  fraction of the total number of sequences in this alignment block. Has been
  superseeded by the \code{gb.*} parameters.}
  \item{\code{block.max.dist}}{Numeric, ranging between 0 and 1.
  \code{block.max.dist} gives the maximum genetic distance (measured as the
  fraction of divergent nucleotide positions) allowed in a sequence alignment
  block. The alignment of individual marker is iteratively broken into
  smaller blocks until this condition is met with.}
  \item{\code{min.n.seq}}{Numeric, the minimum number of sequences required
  for an alignment block. Alignment blocks with less than \code{min.n.seq}
  are dropped from the output.} \item{\code{max.mad}}{Numeric, giving the
  treshold value for the assessment of saturation: alignments with a median
  average distance (MAD) of \code{max.mad} or greater will be broken into
  blocks. The default value has been estimated with simulation by Smith et
  al. (2009).} \item{\code{gb1}}{Parameters for masking of alignment blocks
  with \code{\link{gblocks}}.} \item{\code{gb2}}{Parameters for masking of
  alignment blocks with \code{\link{gblocks}}.} \item{\code{gb3}}{Parameters
  for masking of alignment blocks with \code{\link{gblocks}}.}
  \item{\code{gb4}}{Parameters for masking of alignment blocks with
  \code{\link{gblocks}}.} \item{\code{gb5}}{Parameters for masking of
  alignment blocks with \code{\link{gblocks}}.} }
}

\references{
Smith, S.A., J.M. Beaulieu, and M.J. Donoghue. 2009.
  Mega-phylogeny approach for comparative biology: an alternative to
  supertree and supermatrix approaches. \emph{BMC Evolutionary Biology}
  \bold{9}:37.
}
\seealso{
\code{\link{megapteraProj}} for creating a megaptera project.
}
