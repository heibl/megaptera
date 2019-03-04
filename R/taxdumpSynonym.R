## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2019-03-04)

#' @title Utilities for NCBI Taxdump
#' @description Manage Synonyms
#' @param tax A data frame in parent-child format.
#' @param binomials A vector of mode \code{"character"} giving one or more taxon
#'   names. The first name is considered to be the accepted name.
#' @param keep.acc Logical, only applies when the accepted names in \code{tax}
#'   and \code{binomial} are not identical: \code{TRUE} means the accepted name
#'   from \code{tax} is kept, whereas the accepted name of \code{binomials} is
#'   taken if \code{keep.acc = FALSE}.
#' @param keep.syn Logical, should synonyms present in \code{tax} be kept?
#' @param add.syn Logical, should synonyms not present in \code{tax} be added?
#' @param quiet Logical, if \code{FALSE}, progress is printed on screen for
#'   debugging.
#' @seealso \code{\link{taxdumpAddNode}}, \code{\link{taxdumpChildren}},
#'   \code{\link{taxdumpMRCA}}, \code{\link{taxdumpSubset}},
#'   \code{\link{taxdump2phylo}}
#' @note \code{add.syn} is only partly implemented and not tested, because under
#'   the current use (adding synonyms to NCBI taxonomy) it does not make sense.
#' @export

taxdumpSynonym <- function(tax, binomials, keep.acc = FALSE, 
                           keep.syn = TRUE, add.syn = TRUE, quiet = TRUE){
  
  ## CHECKS
  id <- duplicated(binomials)
  if (any(id)){
    stop("binomial '", binomials[id], "' is duplicated")
  }
  
  ## Assume that first element of 'binomials' is
  ## the accepted name, all others are synonyms
  
  if (!quiet){
    cat("Accepted name:", binomials[1], ".. ")
  }
 
  ## 1. Gather information
  ## a. binomials present in tax
  ## ---------------------
  bb <- tax[tax$taxon %in% binomials , ]
  if (is.null(bb)){
    cat("WARNING: not present (incl. synonyms)\n")
    return(tax)
  }
  status_b <- c("accepted", rep("synonym", length(binomials) - 1))
  status_b <- status_b[binomials %in% tax$taxon]
  status_b <- status_b[match(binomials[binomials %in% tax$taxon], bb$taxon)]
  bb <- cbind(bb, status_b)
  
  ## b. binomials absent from tax
  ## ----------------------------
  id <- !binomials %in% tax$taxon
  if (any(id)){
    temp <- data.frame(
      parent_id = bb$parent_id[1],
      id = bb$id[1],
      taxon = binomials[id],
      rank = "species",
      status = NA,
      status_b = c("accepted", rep("synonyms", length(binomials) - 1))[id],
      stringsAsFactors = FALSE)
    bb <- rbind(bb, temp)
  }
  
  ## c. binomials absent from binomials
  ## ----------------------------------
  temp <- tax[tax$id %in% bb$id, ]
  temp <- temp[!temp$taxon %in% bb$taxon, ]
  if (nrow(temp)){
    temp <- cbind(temp, status_b = NA)
    bb <- rbind(bb, temp)
  }
  
  ## 2. EVALUATION + ADJUSTMENT
  ## a. Only intersection present - nothing to do
  ##    -----------------------------------------
  if (nrow(bb) == 1){
    cat("nothing to do")
    return(tax)
  } 
  
  ## These are the names in tax
  ## tax[tax$taxon %in% bb$taxon, ]
  
  ## b. add missing names (as synonyms, ...)
  ## -----------------------
  id <- is.na(bb$status)
  if (add.syn & any(id)){
    new_taxa <- bb[id, 1:5]
    new_taxa$status <- "synonym"
    tax <- rbind(tax, new_taxa)
  }
  
  ## b. only one name accepted?
  ##    -----------------------
  if (keep.acc){
    stop("implement me!")
  } else {
    id <- grep("accepted", bb$status_b)
    tax$status[tax$taxon %in% bb$taxon[-id]] <- "synonmy"
    tax$status[tax$taxon %in% bb$taxon[id]] <- "scientific name"
  }
  
  ## c. identical parent id?
  ##    --------------------
  if (length(unique(bb$parent_id)) > 1){
    id <- grep("accepted", bb$status_b)
    new_pid <- bb$parent_id[id]
    
    ## Check if the PID's that will be overwritten,
    ## lead to other accepted congenerics. If not, these
    ## parents' status will be set to 'synonym'
    old_pid <- setdiff(bb$parent_id, new_pid)
    for (i in old_pid){
      cg <- tax[tax$parent_id == i & tax$status == "scientific name", ]
      cg <- cg[!cg$taxon %in% bb$taxon, ]
      if (!nrow(cg)){
        tax$status[tax$id == i] <- "synonym"
      }
    }
    tax$parent_id[tax$taxon %in% bb$taxon] <- new_pid
  }
  
  ## d. identical id?
  ##    Note: this can lead to orphaned infraspecific lineages
  ##    ------------------------------------------------------
  if (length(unique(bb$id)) > 1){
    id <- grep("accepted", bb$status_b)
    tax$id[tax$taxon %in% bb$taxon] <- bb$id[id]
  }
  tax
}