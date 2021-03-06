## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2019-02-25)

#' @title Utilities for NCBI Taxdump
#' @description Manage Synonyms
#' @param tax A data frame in parent-child format.
#' @param binomials A vector of mode \code{"character"} giving one or more taxon
#'   names. The first name is considered to be the accepted name, the.
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

taxdumpManageSynonym <- function(tax, binomials, keep.syn = TRUE, add.syn = FALSE, quiet = TRUE){
  
  ## Assume that first element of 'binomials' is
  ## the accepted name, all others are synonyms
  
  if (!quiet){
    cat("Accepted name:", binomials[1], ".. ")
  }
  id <- which(binomials %in% tax$taxon)
  if (!length(id)){
    cat("WARNING: not present (incl. synonyms)\n")
    return(tax)
  }
  
  ## 1. Add accepted name (if it is missing)
  ## ---------------------------------------
  if (min(id) > 1){

    ## Insert missing accepted name
    ## ----------------------------
    parent <- strip.spec(binomials[1])
    if (!parent %in% tax$taxon){
      parent <- sapply(binomials[id], taxdumpHigherRank, x = tax, rank = "family")
      ## DIRTY, DIRTY!
      parent <- table(parent)
      parent <- names(parent)[which.max(parent)]
    } 
    tax <- taxdumpAddNode(taxon = binomials[1], x = tax, rank = "species", parent = parent)
    tax$id[tax$taxon == binomials[1]]
    cat("ADDED AS", parent, "\n")
  } else {
    tax$status[tax$taxon == binomials[1]] <- "scientific name"
    ## Convert rank from subsp. to sp. if necessary
    if (grep("^[[:upper:]][[:lower:]-]+ [[:lower:]-]+$", tax$taxon[tax$taxon == binomials[1]])){
      tax$rank[tax$taxon == binomials[1]] <- "species"
    }
    
    accepted <- tax[tax$taxon == binomials[1], ]
    
    ## Let's check if the name is tied to the right genus name
    ## -------------------------------------------------------
    genus <- taxdumpHigherRank(tax, binomials[1], "genus")
    current_genus <- strip.spec(binomials[1])
    current_parent <- tax$taxon[tax$id == accepted$parent_id & tax$status == "scientific name"]
    
    if (genus != current_genus){
      
      if (nrow(tax[tax$taxon %in% current_genus & tax$rank %in% "genus", ])){
        
        ## Case 1: The corresponding genus name is a synonym and has to be "resurrected"
        if (tax$status[tax$taxon %in% current_genus & tax$rank %in% "genus"] == "scientific name"){
          new_parent_id  <- tax$id[tax$taxon %in% current_genus & tax$rank %in% "genus"]
        } else {
          tax$status[tax$taxon %in% current_genus & tax$rank %in% "genus"] <- "scientific name"
          new_parent_id <- max(tax[, c("id", "parent_id")]) + 1
          tax$id[tax$taxon %in% current_genus & tax$rank %in% "genus"] <- new_parent_id
        }
        tax$parent_id[tax$taxon == binomials[1]] <- new_parent_id
      } else {
        
        ## Case 2: The corresponding genus name must be inserted
        ## a) genus einsetzen
        new_id <- max(tax[, c("id", "parent_id")]) + 1
        new_pid <- tax$parent_id[tax$taxon == genus & tax$rank == "genus"]
        tax <- rbind(tax,
                     c(parent_id = new_pid,
                       id = new_id,
                       taxon = current_genus,
                       rank = "genus",
                       status = "scientific name"))
        ## b) parent_id bei accepted setzen
        tax$parent_id[tax$id == accepted$id] <- new_pid
      }
    }
    if (current_parent != current_genus){
      cat("I am alive!")
      ## The corresponding genus name is a synonym and has to be "resurrected"
      new_parent_id  <- tax$id[tax$taxon %in% current_genus & tax$rank %in% "genus"]
      tax$parent_id[tax$taxon == accepted$taxon] <- new_parent_id 
    }
    cat("OK\n")
  }
  
  ## 2. Set ID/PARENT_ID of synonyms to ID/PARENT_ID of accepted name
  ## ----------------------------------------------------------------
  tax_id <- tax$id[tax$taxon == binomials[1]]
  tax$id[tax$taxon %in% binomials[id]] <- tax_id
  # tax$parent_id[tax$id == tax_id] <- tax$id[tax$taxon == strip.spec(binomials[1])]
  # The line above only works when there are now homonymic subgenera (e.g., 
  # Equus = genus and subgenus), in this case, do:
  tax$parent_id[tax$id == tax_id] <- tax$parent_id[tax$taxon == binomials[1]]
  
  
  
  ## 3. Make sure that all other names (also those coming from 'tax')
  ##    will be marked 'synonym'
  ## ---------------------------
  syn <- tax$taxon[tax$id == tax_id]
  syn <- setdiff(syn, binomials[1])
  tax$status[tax$taxon %in% syn] <- "synonym"
  
  ## 4. Manage synonyms
  ## ------------------
  if (keep.syn){
    if (add.syn){
     syn_to_add <- binomials[-id]
     if (length(syn_to_add)){
       family <- sapply(binomials[1], taxdumpHigherRank, x = tax, rank = "family")
       for (i in syn_to_add){
         tax <- taxdumpAddNode(taxon = i, x = tax, rank = "species", parent = family)
         ## CAUTION: Will set status to "scientific name"
       }
     }
    } else {
      
    }
  } else {
    
    ## Delete synonyms (if there are any and 'keep.syn' == FALSE)
    if (length(id[id > 1])){
      tax <- tax[!tax$taxon %in% binomials[id[id > 1]], ]
    }
  }
  tax
}