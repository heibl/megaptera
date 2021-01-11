## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2019-11-05)


#' @title Convert SQL Query Result to DNAbin
#' @description Convert a data frame as returned by \code{dbGetQuery} to an
#'   object of class \code{\link{DNAbin}}, possibly with confidence scores.
#' @param pg A data frame as returned by \code{dbGetQuery}.
#' @param label A vector of mode \code{"character"} giving the attributes for
#'   constructing the sequences labels, e.g. \code{"taxon"}, \code{"acc"}, or both.
#' @param confid.scores A vector of mode \code{"character"} specifying if and
#'   how confidence scores will be returned. Use \code{"ignore"} to omit
#'   confidence scores. Other options are \code{"all"}, \code{"row.means"}, and
#'   \code{"col.means"}. Confidence scores are returned as an attribute and can
#'   be accessed by \code{attr(obj, "cs")}.
#' @param row.confid A real number in the interval \code{[0, 1]} giving
#'   the confidence threshold for alignment rows (i.e. taxa, sequences); only rows  scoring equal
#'   or greater to \code{min.confid} will be selected.
#' @param col.confid A real number in the interval \code{[0, 1]} giving
#'   the confidence threshold for alignment columns (i.e. nucletide positions); only rows
#'   (i.e. taxa, sequences) or columns  scoring equal
#'   or greater to \code{min.confid} will be selected.
#' @importFrom ape as.DNAbin
#' @export

pg2DNAbin <- function(pg, label, confid.scores, col.confid = 0, row.confid = 0){
  
  ## CHECKS
  ## ------
  if (!is.data.frame(pg)) stop("'pg' must be of class 'data.frame'")
  id <- c("taxon", "sequence", "reliability") %in% names(pg)
  if (!all(id)) stop("'pg' is malformatted")
  confid.scores <- match.arg(confid.scores, c("ignore", "all", "row.means", "col.means")) 
  
  ## Parse DNA sequences
  ## -------------------
  obj <- strsplit(pg$sequence, split = "")
  names(obj) <- apply(pg[, label, drop = FALSE], 1, paste, collapse = "_")
  names(obj) <- gsub(" ", "_", names(obj))
  obj <- as.DNAbin(obj)
  if (length(unique(sapply(obj, length))) == 1){
    obj <- as.matrix(obj)
  } else {
    # ## If sequences are not aligned, they cannot
    # ## have reliability scores assigned to them
    # return(obj)
  }
  
  ## If reliability scores are not available,
  ## stop here and return sequences
  ## ------------------------------
  test <- sapply(pg$reliability, nchar) == 1
  if (any(is.na(test)) | confid.scores == "ignore"){
    return(obj) ## reliability is NA -> not yet aligned
  }
  if (any(test) | confid.scores == "ignore"){
    return(obj) ## reliability is a single number -> only alignment
  }
  
  ## Parse reliability score. Note this always returns a matrix, 
  ## even if only row or column scores are available
  ## -----------------------------------------------
  rel_score <- strsplit(pg$reliability, split = " ")
  rel_score <- lapply(rel_score, as.numeric)
  rel_score <- do.call(rbind, rel_score)
  if (!identical(dim(obj), dim(rel_score))){
    stop("mismatch between alignment and reliability scores")
  }
  
  ## Select rows (= sequences) >= row.confid
  ## ---------------------------------------
  if (row.confid > 0){
    id <- rowMeans(rel_score) >= row.confid
    obj <- obj[id, ]
    rel_score <- rel_score[id, ]
    if (!nrow(obj)) stop("please use a smaller value for 'row.confid'")
  }
  
  ## Select columns (= nucleotide positions) >= col.confid
  ## -----------------------------------------------------
  if (col.confid > 0){
    id <- colMeans(rel_score) >= col.confid
    obj <- obj[, id]
    rel_score <- rel_score[, id]
    if (!ncol(obj)) stop("please use a smaller value for 'col.confid'")
  }
  
  ## Simplify condidence scores if desired
  ## -------------------------------------
  if (confid.scores == "row.means"){
    rel_score <- rowMeans(rel_score, na.rm = TRUE)
  }
  if (confid.scores == "col.means"){
    rel_score <- colMeans(rel_score, na.rm = TRUE)
  }
  
  attr(obj, "cs") <- rel_score
  obj
}
