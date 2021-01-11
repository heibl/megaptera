## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-09-17)

#' @title Class "megapteraPars"
#' @description S4 Class for pipeline parameters of \code{\link{megapteraProj}}.
#' @slot data.path A character string giving the path to the directory there all
#'   data and results will be stored (see \code{\link{megapteraInit}}).
#' @slot gb.seq.download A character string defining how sequences should be
#'   downloaded from GenBank Nucleotide; can be \code{"eutils"} or \code{"ftp"}.
#' @slot debug.level Numeric, a number between 0 and 5, determining the
#'   pipeline's verbosity. \emph{This parameter under development and should not
#'   be changed by the user.}
#' @slot parallel Logical: if \code{TRUE}, several steps in the pipeline will be
#'   run in parallel, otherwise all steps are serial.
#' @slot cpus Numerical: if \code{TRUE}, several steps in the pipeline will be
#'   run in parallel, otherwise all steps are serial.
#' @slot cluster.type A character string: if \code{TRUE}, several steps in the
#'   pipeline will be run in parallel, otherwise all steps are serial.
#' @slot update.seqs \emph{Currently unused}.
#' @slot retmax Numeric, giving the batch size when downloading sequences from
#'   the Entrez History server (default: 500).
#' @slot max.gi.per.spec Numeric, giving the maximum number of sequences that
#'   will be used per species. Can be used to avoid model organism (e.g., rice,
#'   \emph{Drosophila}, ...) cluttering up the pipeline with thousands of
#'   sequences (default: 10).
#' @slot max.bp Numeric, the maximal length of DNA sequences in base pairs to be
#'   included in the alignment. The upper limit is determined by the alignment
#'   program and the specific alignment and can only be determined by
#'   trial-and-error (default: 5000).
#' @slot reference.max.dist \emph{Currently unused}.
#' @slot min.seqs.reference \emph{Currently unused}.
#' @slot fract.miss Numeric, ranging between 0 and 1. To avoid long stretches of
#'   only a few sequences at the beginning and the ending of an alignment block
#'   a minimum required number of sequences can be set as a fraction of the
#'   total number of sequences in this alignment block. Has been superseeded by
#'   the \code{gb.*} parameters.
#' @slot filter1 \emph{Currently unused}.
#' @slot filter2 \emph{Currently unused}.
#' @slot filter3 \emph{Currently unused}.
#' @slot filter4 \emph{Currently unused}.
#' @slot block.max.dist Numeric, ranging between 0 and 1. \code{block.max.dist}
#'   gives the maximum genetic distance (measured as the fraction of divergent
#'   nucleotide positions) allowed in a sequence alignment block. The alignment
#'   of individual marker is iteratively broken into smaller blocks until this
#'   condition is met with.
#' @slot min.n.seq Numeric, the minimum number of sequences required for an
#'   alignment block. Alignment blocks with less than \code{min.n.seq} are
#'   dropped from the output.
#' @slot max.mad Numeric, giving the treshold value for the assessment of
#'   saturation: alignments with a median average distance (MAD) of
#'   \code{max.mad} or greater will be broken into blocks. The default value has
#'   been estimated with simulation by Smith et al. (2009).
#' @slot gb1 Parameters for masking of alignment blocks with
#'   \code{\link{gblocks}}.
#' @slot gb2 Parameters for masking of alignment blocks with
#'   \code{\link{gblocks}}.
#' @slot gb3 Parameters for masking of alignment blocks with
#'   \code{\link{gblocks}}.
#' @slot gb4 Parameters for masking of alignment blocks with
#'   \code{\link{gblocks}}.
#' @slot gb5 Parameters for masking of alignment blocks with
#'   \code{\link{gblocks}}. #' @details The pipeline's verbosity can be
#'   fine-tuned with \code{debug.level}: \tabular{ll}{ 0 \tab No progess and
#'   diagnostic messages\cr 1 \tab Messages on screen\cr 2 \tab Messages logged
#'   to file\cr 3 \tab Messages on screen and logged to file\cr 4 \tab Same as
#'   3, in addition current data is saved as .rda object in case of a
#'   foreseeable error\cr 5 \tab Same as 4, in addition current data is always
#'   saved\cr }
#' @references  Smith, S.A., J.M. Beaulieu, and M.J. Donoghue. 2009.
#'   Mega-phylogeny approach for comparative biology: an alternative to
#'   supertree and supermatrix approaches. \emph{BMC Evolutionary Biology}
#'   \bold{9}:37.
#' @seealso \code{\link{megapteraProj}} for creating a megaptera project.

setClass("megapteraPars", 
         representation = list(
           data.path = "character", 
           gb.seq.download = "character", # stepB
           debug.level = "numeric", # all steps
           parallel = "logical",
           cpus = "numeric",
           cluster.type = "character",
           update.seqs = "character", # step B
           retmax = "numeric", # step B
           max.gi.per.spec = "numeric", # step B
           max.bp = "numeric", # step B
           reference.max.dist = "numeric", # step C
           min.seqs.reference = "numeric", # step D
           fract.miss = "numeric", # step G
           filter1 = "numeric",  # step G
           filter2 = "numeric", # step G
           filter3 = "numeric", # step G
           filter4 = "numeric",  # step G
           block.max.dist = "numeric", # step G
           min.n.seq = "numeric", # step G
           max.mad = "numeric", # step H
           gb1 = "numeric", # step F + G
           gb2 = "numeric", # step F + G
           gb3 = "numeric", # step F + G
           gb4 = "numeric", # step F + G
           gb5 = "character") # step F + G
)
