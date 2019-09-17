## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-09-11)

#' @title Step B: Search and Download Sequences
#' @description For any given project (see \code{\link{megapteraProj}}), 
#' this step searches the Nucleotide database on GenBank, downloads all 
#' sequences and stores them in a postgreSQL database table.
#' @details All accessions are stored under their species name as appearing 
#' in the \emph{organism} field at GenBank, but information about infrageneric 
#' ranks is stripped off the taxon names before they are stored in the database.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param update.seqs A character string determining the behaviour of 
#' \code{stepB} when it is run repeatedly on the same locus/taxon 
#' combination: \code{"no"} means that only sequences that have been made 
#' available on GenBank since the last execution \code{stepB} will be 
#' downloaded; \code{"all"} means the database table with all sequences 
#' will be removed and all sequences will be downloaded de novo. Naturely, 
#' this option has no effect if \code{stepB} is called for the first time.
#' @return None. \code{stepB} is called for its side effects: (1) strings of 
#' DNA sequences with attribute data are stored in a pgSQL database; (2) a log 
#' file is written to the current working directory.
#' @references \url{http://www.ncbi.nlm.nih.gov/books/NBK25501/}
#' @seealso \code{\link{megapteraProj}}; \code{\link{stepA}} for the preceeding 
#' and \code{\link{stepC}} for the subsequent step; \code{\link{stepBX}} for the 
#' addition of external sequences to the database.
#' @export
#' @import DBI

stepB <- function(x, update.seqs = "no"){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("undefined locus not allowed")
  if (!url.exists("https://eutils.ncbi.nlm.nih.gov"))
    stop("internet connection required for stepB")
  
  ## Switch between E-utilities or FTP to download
  ## sequences from GenBank::Nucleotide
  ## ----------------------------------
  if (x@params@gb.seq.download == "ftp"){
    stepB_ftp(x = x, update.seqs = update.seqs)
  } else {
    stepB_eutils(x = x, update.seqs = update.seqs)
  }
}
