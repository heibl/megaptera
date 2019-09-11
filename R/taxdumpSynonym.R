## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2019-07-18)

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
  
  ## Columns 'parent_id' and 'id' must be numeric
  ## -------------------------------------------
  tax$parent_id <- as.numeric((tax$parent_id))
  tax$id <- as.numeric((tax$id))
  
  ## Assume that first element of 'binomials' is
  ## the accepted name, all others are synonyms
  
  if (!quiet){
    cat("Accepted name:", binomials[1], ".. ")
  }
 
  #################################################
  ## PART A. Gather information
  #################################################
  
  ## A1. binomials present in tax
  ## ----------------------------
  bb <- tax[tax$taxon %in% binomials , ]
  if (is.null(bb)){
    cat("WARNING: not present (incl. synonyms)\n")
    return(tax)
  }
  status_b <- c("accepted", rep("synonym", length(binomials) - 1))
  status_b <- status_b[binomials %in% tax$taxon]
  status_b <- status_b[match(binomials[binomials %in% tax$taxon], bb$taxon)]
  bb <- cbind(bb, status_b)
  
  ## A2. binomials absent from tax
  ## -----------------------------
  id <- !binomials %in% tax$taxon
  if (any(id)){
    temp <- data.frame(
      parent_id = bb$parent_id[1],
      id = bb$id[1],
      taxon = binomials[id],
      rank = "species",
      status = NA,
      status_b = c("accepted", rep("synonym", length(binomials) - 1))[id],
      stringsAsFactors = FALSE)
    bb <- rbind(bb, temp)
  }
  
  ## A3. binomials absent from binomials
  ## -----------------------------------
  temp <- tax[tax$id %in% bb$id, ]
  temp <- temp[!temp$taxon %in% bb$taxon, ]
  if (nrow(temp)){
    temp <- cbind(temp, status_b = NA)
    bb <- rbind(bb, temp)
  }
  bb_archive <- bb
  
  #################################################
  ## PART B. EVALUATION + ADJUSTMENT
  #################################################
  
  ## B1. Only intersection present - nothing to do
  ## ---------------------------------------------
  if (nrow(bb) == 1){
    cat("nothing to do")
    return(tax)
  } 
  
  ## These are the names in tax
  ## tax[tax$taxon %in% bb$taxon, ]
  
  ## B2. add missing names (as synonyms; if necessary 
  ##     their status will be modified later)
  ## -------------------------------------------------
  id <- is.na(bb$status)
  if (add.syn & any(id)){
    new_taxa <- bb[id, 1:5]
    new_taxa$status <- "synonym"
    tax <- rbind(tax, new_taxa)
  }
  
  ## B3. 'binomials' can match with more than one accepted
  ##     name in 'tax'; select the one that is accepted 
  ##     according to 'binomials' or - if not possible -
  ##     select one randomly
  ## -----------------------------------------------------
  accepted_t <- which(bb$status == "scientific name")
  accepted_b <- which(bb$status_b == "accepted")
  if (length(accepted_t) > 1){
    if (length(intersect(accepted_t, accepted_b))) {
      id <- bb$id[setdiff(accepted_t, accepted_b)]
    } else {
      id <- bb$id[accepted_t][1]
    }
    tax$status[tax$id == id] <- "synonym"
    bb$status[bb$id == id] <- "synonym"
  }
  
  ## B4. only one name accepted!
  ## ----------------------------
  if (keep.acc){
    stop("implement me!")
  } else {
    id <- grep("accepted", bb$status_b)
    tax$status[tax$taxon %in% bb$taxon[-id]] <- "synonym"
    tax$status[tax$taxon %in% bb$taxon[id]] <- "scientific name"
    
    ## Make sure that for recombination the new genus name is
    ## available and accepted
    ## ----------------------
    genus_t <- unique(strip.spec(bb$taxon[bb$status %in% "scientific name"]))
    genus_b <- strip.spec(bb$taxon[bb$status_b %in% "accepted"])
    
    if (length(c(genus_t, genus_b)) > 2) stop("debug me (code tb)")
    
    ## This is what happens in case of a recombination:
    if (genus_t != genus_b){
      ## Is genus available?
      id <- genus_b %in% tax$taxon
      if (length(which(id)) > 2) stop("debug me!")
      if (id){
        ## Is 'genus_b' a synonym?
        id <- tax$taxon %in% genus_b
        if (tax$status[id] == "synonym"){
          ## Set status to "accepted"
          tax$status[id] <- "scientific name"
          ## Assign new ID 
          nid <- max(tax[, c("id", "parent_id")]) + 1
          tax$id[id] <- nid
          tax$parent_id[tax$taxon %in% bb$taxon] <- nid
          bb$parent_id[] <- nid
        } else {
          ## 'genus_b' is available and accepted, 
          ## so use its ID
          tax$parent_id[tax$taxon %in% bb$taxon] <- tax$id[id]
          bb$parent_id[] <- tax$id[id]
        }
      } else {
        ## genus_b is not available
        ## create its entry and link to children
        ## ASSUPMTION: 'genus_t' and 'genus_b' share same parent
        nid <- max(tax[, c("parent_id", "id")]) + 1
        new_genus <- data.frame(
          parent_id = tax$parent_id[tax$taxon == genus_t], 
          id = nid,
          taxon = genus_b,
          rank = "genus",
          status = "scientific name",
          stringsAsFactors = FALSE)
        tax$parent_id[tax$taxon %in% bb$taxon] <- nid 
        bb$parent_id[] <- nid
        tax <- rbind(tax, new_genus)
        
      }
      ## Check if 'genus_t' results childless
      test <- taxdumpChildren(tax, genus_t, query.rank = "genus")
      test <- test[test$rank != "genus", ]
      if (!nrow(test)){
        ## ... yes it is, so consider it a synonym of 'genus_b'
        opid <- tax$id[tax$taxon == genus_t & tax$rank == "genus"]
        npid <- tax$id[tax$taxon == genus_b & tax$rank == "genus"]
        tax[tax$id == opid, c("id", "status")] <- data.frame(npid, "synonym",
                                                             stringsAsFactors = FALSE)
      }
    }
  }
  
  ## B5. identical parent id?
  ## ----------------------------------
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
  
  ## B6. identical id?
  ## Note: this can lead to orphaned infraspecific lineages
  ## ------------------------------------------------------
  if (length(unique(bb$id)) > 1){
    id <- grep("accepted", bb$status_b)
    tax$id[tax$taxon %in% bb$taxon] <- bb$id[id]
  }
  
  if (any(is.na(tax$parent_id)) | any(is.na(tax$id)))
    stop("debug me! (NAs generated")
  
  tax
}