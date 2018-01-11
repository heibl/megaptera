## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-01-11)

#' @title Create an Object of Class "megapteraPars"
#' @description S4 Class for parameters of a megaptera project pipeline, as stored in \code{\link{megapteraProj}}.
#' @param ... Arguments in \code{tag = value} form. The tags must come from the
#'   names of the parameters described in the `Pipeline Parameters' section.
#' @section Pipeline Parameters: 
#' \describe{ 
#' \item{\code{debug.level}}{Numeric,
#'   a number between 0 and 5, determining the pipeline's verbosity. \emph{This
#'   parameter is under development and should not be changed by the user.}} 
#' \item{\code{parallel}}{Logical: if \code{TRUE}, several steps in the pipeline will
#'   be run in parallel, otherwise all steps are serial.}
#' \item{\code{cpus}}{Numerical: if \code{TRUE}, several steps in the pipeline will be
#'   run in parallel, otherwise all steps are serial.}
#' \item{\code{cluster.type}}{A character string: if \code{TRUE}, several steps in the
#'   pipeline will be run in parallel, otherwise all steps are serial.}
#' \item{\code{update.seqs}}{\emph{Currently unused}.}
#' \item{\code{retmax}}{Numeric, giving the batch size when downloading sequences from
#'   the Entrez History server (default: 500).}
#' \item{\code{max.gi.per.spec}}{Numeric, giving the maximum number of sequences that
#'   will be used per species. Can be used to avoid model organism (e.g., rice,
#'   \emph{Drosophila}, ...) cluttering up the pipeline with thousands of
#'   sequences (default: 10).}
#' \item{\code{max.bp}}{Numeric, the maximal length of DNA sequences in base pairs to
#'   be included in the alignment. The upper limit is determined by the
#'   alignment program and the specific alignment and can only be determined by
#'   trial-and-error (default: 5000).}
#' \item{\code{reference.max.dist}}{\emph{Currently unused}.}
#' \item{\code{min.seqs.reference}}{\emph{Currently unused}.}
#' \item{\code{fract.miss}}{Numeric, ranging between 0 and 1. To avoid long stretches
#'   of only a few sequences at the beginning and the ending of an alignment
#'   block a minimum required number of sequences can be set as a fraction of
#'   the total number of sequences in this alignment block. Has been superseeded
#'   by the \code{gb.*} parameters.}
#' \item{\code{block.max.dist}}{Numeric, ranging between 0 and 1. \code{block.max.dist}
#'   gives the maximum genetic distance (measured as the fraction of divergent
#'   nucleotide positions) allowed in a sequence alignment block. The alignment
#'   of individual marker is iteratively broken into smaller blocks until this
#'   condition is met with.}
#' \item{\code{min.n.seq}}{Numeric, the minimum number of sequences required for an
#'   alignment block. Alignment blocks with less than \code{min.n.seq} are
#'   dropped from the output.}
#' \item{\code{max.mad}}{Numeric, giving the treshold value for the assessment of
#'   saturation: alignments with a median average distance (MAD) of
#'   \code{max.mad} or greater will be broken into blocks. The default value has
#'   been estimated with simulation by Smith et al. (2009).}
#' \item{\code{gb1}}{Parameters for masking of alignment blocks with
#'   \code{\link{gblocks}}.}
#' \item{\code{gb2}}{Parameters for masking of alignment blocks with
#'   \code{\link{gblocks}}.}
#' \item{\code{gb3}}{Parameters for masking of alignment blocks with
#'   \code{\link{gblocks}}.}
#' \item{\code{gb4}}{Parameters for masking of alignment blocks with
#'   \code{\link{gblocks}}.}
#' \item{\code{gb5}}{Parameters for masking of alignment blocks with
#'   \code{\link{gblocks}}.}
#' }
#' @references  Smith, S.A., J.M. Beaulieu, and M.J. Donoghue. 2009.
#'   Mega-phylogeny approach for comparative biology: an alternative to
#'   supertree and supermatrix approaches. \emph{BMC Evolutionary Biology}
#'   \bold{9}:37.
#' @details The pipeline's verbosity can be fine-tuned with \code{debug.level}:
#'   \tabular{ll}{ 0 \tab No progess and diagnostic messages\cr 1 \tab Messages
#'   on screen\cr 2 \tab Messages logged to file\cr 3 \tab Messages on screen
#'   and logged to file\cr 4 \tab Same as 3, in addition current data is saved
#'   as .rda object in case of a foreseeable error\cr 5 \tab Same as 4, in
#'   addition current data is always saved\cr }
#' @seealso \code{\link{megapteraProj}} for creating a megaptera project.
#' @examples
#' megapteraPars()
#' @include megapteraPars-class.R
#' @importFrom methods new slotNames
#' @export

"megapteraPars" <- function(...){
  
  params <- list(debug.level = 1,
                 parallel = FALSE,
                 cpus = 0,
                 cluster.type = "none",
                 update.seqs = "all", 
                 retmax = 500,
                 max.gi.per.spec = 10, 
                 max.bp = 5000, 
                 reference.max.dist = 0.25,
                 min.seqs.reference = 10,
                 fract.miss = .25, 
                 block.max.dist = .5, # step G
                 min.n.seq = 5, # step G
                 max.mad = .01,
                 gb1 = .5, # step F + G
                 gb2 = .5, # step F + G
                 gb3 = 9999, # step F + G
                 gb4 = 2, # step F + G
                 gb5 = "a" # step F + G
  )
  
  args <- list(...)
  notDef <- setdiff(names(args), names(params))
  if ( length(notDef) ) 
    stop ("parameter '", notDef[1], "' is not defined", sep = "")
  
  id <- match(names(args), names(params))
  params[id] <- args
  
  if (params$parallel & params$cpus == 0){
    stop("number of CPUs must be given")
  }
  if (params$parallel & params$cluster.type == "none"){
    stop("type of cluster must be given: 'SOCK', 'MPI', 'PVM', or 'NWS'")
  }
  
  new("megapteraPars",
      debug.level = params$debug.level,
      parallel = params$parallel,
      cpus = params$cpus,
      cluster.type = params$cluster.type,
      update.seqs = params$update.seqs, 
      retmax = params$retmax,
      max.gi.per.spec = params$max.gi.per.spec, 
      max.bp = params$max.bp,
      reference.max.dist = params$reference.max.dist,
      min.seqs.reference = params$min.seqs.reference,
      fract.miss = params$fract.miss, 
      block.max.dist = params$block.max.dist,
      min.n.seq = params$min.n.seq,
      max.mad = params$max.mad,
      gb1 = params$gb1,
      gb2 = params$gb2,
      gb3 = params$gb3,
      gb4 = params$gb4,
      gb5 = params$gb5
  )
}

setMethod("show",
          signature(object = "megapteraPars"),
          function (object) 
          {
            out <- sapply(slotNames(object), slot, object = object)
            names(out) <- format(names(out), justify = "right")
            out <- paste("\n ", names(out), "=", out)
            cat("MEGAPTERA pipeline parameters:", out)
          }
)