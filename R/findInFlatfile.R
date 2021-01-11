
#' @title Find Taxa and Accession Numbers in GenBank Flatfiles
#' @description Checks if any entries of a taxon are present in a given GenBank
#'   flatfile based on regular expressions. The motivation is to save computation
#'   time by not parsing files that do not include the taxon of interest.
#' @param file A character string giving the name of the GenBank flatfile.
#' @param taxon A character string giving the name of the taxon.
#' @param acc A character string giving the accession number.
#' @return Logical, returns \code{TRUE}, if \code{taxon} is included with one or
#'   more entries in \code{file}.
#' @keywords internal
#' @export

findInFlatfile <- function(file, taxon = NULL, acc = NULL){
  
  ## REMEMBER:
  ## 1. zgrep calls grep directly
  ## 2. Infix operator and '+' must be escaped in normal grep
  ## 3. flag -m 1 means return first match in any given file
  call <- NULL
  if (!is.null(taxon)){
    call <- c(call, paste0("ORGANISM[[:space:]]\\+", taxon))
  }
  if (!is.null(acc)){
    call <- c(call, paste0("ACCESSION[[:space:]]\\+", acc))
  }
  call <- paste(call, collapse = "\\|")
  call <- paste0("'", call, "'")
  call <- paste("zgrep -m 1", call, file)
  suppressWarnings(out <- system(call, intern = TRUE))
  length(out) == 1
}
