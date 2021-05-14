## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2021-05-11)

#' @title Sanity Check for Parent-Child Taxonomic Tables
#' @description Does several sanity checks for taxonomic tables in parent-child format.
#' @param tax A data frame in parent-child format.
#' @param quiet Logical, indicating id diagnostic output to screen should be supressed.
#' @details \code{taxdumpSanity} checks, if
#' \enumerate{
#' \item accepted taxa and their synonyms share the same parent
#' \item accepted taxa are not linked to parent taxa of the same rank
#' \item accepted species are linked to the correct genus name
#' }
#' @return Logical
#' @export

taxdumpSanity <- function(tax, quiet = FALSE){
  
  is_sane <- TRUE
  
  ## Check input format
  if (!inherits(tax, "data.frame")){
    if (!quiet) cat("FATAL: 'tax' is not of class 'data.frame'\n")
    return(FALSE)
  }
  
  
  id <- duplicated(tax)
  if (any(id) ){
    if (!quiet) cat("\nFATAL: 'tax' contains duplicate rows")
    is_sane <- FALSE
  }
  
  if (!"status" %in% names(tax)){
    tax <- data.frame(tax, status = "scientific name",
                      stringsAsFactors = FALSE)
  }
    
  ## Summary
  ## -----
  if (!quiet) cat("Number of taxon names:", nrow(tax))
  id_set <- sort(unique(tax$id))
  if (!quiet) cat("\nNumber of taxon concepts:", length(id_set))
  
  ## Synonyms are defined by having the same id as accepted taxa
  with_syn <- length(unique(tax$id[tax$status == "synonym"]))
  if (!quiet) cat("\nNumber of taxon concepts with synonyms:", with_syn)
  
  if (!quiet) cat("\n\nMake sure that ...")
  
  ## 1. Do synonyms have same parent as accepted taxa?
  ## -------------------------------------------------
  if (!quiet) cat("\n... accepted taxa and their synonyms share the same parent ...")
  parents <- tapply(tax$parent_id, tax$id, function(z) length(unique(z)))
  parents <- names(parents)[parents > 1]
  if (length(parents)){
    cat("\nFATAL:", length(parents), "taxon", 
        ifelse(length(parents) == 1, "concept has", "concepts have"), 
        "more than one parent:", formatSpecList(parents, n.element = 6))
    is_sane <- FALSE
  } else {
    if (!quiet) cat(" OK")
  }
  
  ## 2. Find taxon names that are linked to a name of the same rank
  ## --------------------------------------------------------------
  if (!quiet) cat("\n... accepted taxa are not linked to parent taxa of the same rank ...")
  test <- tax[tax$status == "scientific name", ]
  test$parent_rank <- test$rank[match(test$parent_id, test$id)] 
  ## Find and exclude root
  root <- which(is.na(test$parent_rank))
  if (length(root) > 1) stop ("debug me!")
  test <- test[-root, ] ## This is the root
  test <- test[test$rank != "no rank" & test$parent_rank != "no rank", ]
  id <- test$rank == test$parent_rank
  if (any(id)){
    n <- nrow(test[id, ])
    p <- test[match(test$parent_id[id], test$id), c("taxon", "rank")]
    p <- data.frame(test$taxon[id], p, stringsAsFactors = FALSE)
    names(p)  <- c("child", "parent", "shared rank")
    cat("\nFATAL:", n, "accepted taxon", ifelse(n == 1, "name is", "names are"), 
        "linked to a taxon name of the same rank:", formatDF(p))

    return(FALSE)
  } else {
    if (!quiet) cat(" OK")
  }
  
  ## 3. Find accepted species linked to the wrong genus. This may happen,
  ##    when recombination was done
  ## ------------------------------
  if (!quiet) cat("\n... accepted species are linked to the correct genus name ...")
  accepted <- tax[tax$status == "scientific name", ]
  children <- accepted[accepted$rank == "species", c("parent_id", "id", "taxon")]
  children <- children[-grep(indet.strings(collapse = TRUE), children$taxon), ]
  parents <- accepted[match(children$parent_id, accepted$id), c("id", "taxon", "rank", "status")]
  test <- cbind(parents, children[, -1])
  names(test) <- c("parent_id", "parent_taxon", "parent_rank", "parent_status","child_id", "child_taxon")
  genus <- test[test$parent_rank == "genus", c("parent_taxon", "child_taxon")]
  genus$genus <- strip.spec(genus$child_taxon, mode = "regex")
  id <- genus$parent_taxon != genus$genus
  # head(genus[id, ])
  if (any(id)){
    n <- length(genus$child_taxon[id])
    genus <- genus[id, c("parent_taxon", "child_taxon")]
    names(genus) <- c("genus", "species")
    cat("\nFATAL:", n, "accepted species", ifelse(n == 1, "name is", "names are"), 
        "linked to a non-corresponding genus name:", 
        formatDF(genus))
    is_sane <- FALSE
  } else {
    if (!quiet) cat(" OK")
  }
  test <- test[test$parent_rank != "genus", ]
  if (nrow(test) & !quiet){
    cat("\nNOTE:", nrow(test), 
        "accepted species are not direktly linked to a genus,", 
        "but to one of these ranks", 
        # formatSpecList(sort(unique(test))))
        formatSpecList(unique(test)))
  }
  remove(test)
  
  ## UNDER DEVELOPMENT: MUST BE DONE FOR KINGDOMS SEPARATELY.
  ## 4. There must not be duplicated accepted names of the same rank
  ## ---------------------------------------------------------------
  # cat("\n... there are no duplicated accepted names of same rank ...")
  # metazoa <- taxdumpSubset(tax, "Metazoa")
  # metazoa <- metazoa[metazoa$status == "scientific name" & metazoa$rank != "no rank",
  #             c("taxon", "rank")]
  # id <- duplicated(metazoa)
  # if (any(id) {
  #   metazoa[test, ]
  # } else {
  #   cat(" OK")
  # }
  if (!quiet) cat("\n")
  is_sane
}
