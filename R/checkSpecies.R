## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-10-19)

#' @title Get Information about Species
#' @description Returns information about the 'fate' of a single species along the pipeline.
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @param spec A character string giving the name of a species.
#' @return A data frame, as a side effect a summary message is printed on the screen.
#' @seealso \code{\link{checkMissingSpec}} to create a list of species that were missed/lost during consecutive steps of the pipeline.
#' @import DBI
#' @export

checkSpecies <- function(megProj, spec){
  
  spec <- gsub("_", " ", spec)
  conn <- dbconnect(megProj)
  
  ## A: TAXONOMY
  ## -----------
  tax <- dbReadTaxonomy(megProj)
  
  ## If not found, is it a synonym?
  ## ------------------------------
  # if (!nrow(tax)){
  #   tax <- paste("SELECT spec FROM  taxonomy", 
  #                "WHERE", wrapSQL(spec, "synonym"))
  #   tax <- dbGetQuery(conn, tax)$spec
  #   
  #   if (!length(tax)){
  #     cat("\n'", spec, "' not present in 'taxonomy' table\n", sep = "")
  #     dbDisconnect(conn)
  #     return()
  #   } else {
  #     cat("NOTE:", spec, "is coded as synonym of", tax, "\n")
  #     spec <- tax
  #     tax <- paste <- paste("SELECT * FROM  taxonomy", 
  #                           "WHERE", wrapSQL(tax))
  #     tax <- dbGetQuery(conn, tax)
  #   }
  # }
  tax <- taxdumpLineage(tax, spec)
  if (is.null(tax)){
    cat("\n'", spec, "' not present in 'taxonomy' table\n", sep = "")
  } else {
    tax <- tax[, c("rank", "taxon")]
    tax <- format(tax, justify = "right")
    tax <- paste(tax$rank, tax$taxon, sep = " : ")
    cat("***** TAXONOMY *****", paste0("\n", tax))
  }
  
  
  ## B: SEQUENCES
  ## ------------
  cat("\n\n***** SEQUENCES *****")
  total <- paste("SELECT locus , taxon FROM sequence WHERE", 
                 wrapSQL(spec, "taxon", "="))
  total <- dbGetQuery(conn, total)
  if (!nrow(total)){
    cat("\nNo sequences of species '", spec, 
        "' in database\n", sep = "")
    dbDisconnect(conn)
    return()
  } else {
    cat("\nNumber of sequences in database:", nrow(total))
  }
  
  ## C: LOCI
  ## -------
  if (any(!is.na(total$locus))){
    nass <- table(total$locus) ## number assigned
    nloc <- length(nass)
    nloc <- ifelse(nloc == 1, paste(nloc, "locus"), paste(nloc, "loci"))
    nass <- sum(nass)
    nass <- ifelse(nass == 1, paste(nass, "of these has"), paste(nass, "of these have"))
    cat("\n", nass ," been assigned (by BLAST) to ", nloc,  "\n", sep = "")
  } else {
    cat("\nNone of these has been assigned (by BLAST) to a particular locus")
    dbDisconnect(conn)
    return()
  }
  
  tab <- dbReadLocus(megProj)
  tab <- tab[rownames(tab) == spec, ]

  loci <- colnames(tab)[tab > 0]
  loci <- sort(unique(gsub("(gb|sel)_", "", loci)))
  loci <- gsub("(^[[:digit:]])", "_\\1", loci) ## _5_8S
  
  ACC <- data.frame()
  
  
  # tabb <- tab[, grep(loci[i], colnames(tab)), drop = FALSE]
  SQL <- paste("SELECT locus, acc, length, evalue, qcovs",
               "FROM sequence",
               "WHERE", wrapSQL(spec, term = "taxon"))
  acc <- dbGetQuery(conn, SQL)
  acc <- acc[order(acc$locus, acc$evalue, acc$qcovs, 
                   decreasing = c(FALSE, FALSE, TRUE),
                   method = "radix"), ]
  rownames(acc) <- NULL
  acc <- split(acc, f = acc$locus)
  
  ## screen output
  ## -------------
  for (i in loci){
    cat("\n*** LOCUS : ", i, " ***", sep = "")
    lc <- acc[names(acc) == i][[1]]
    cat("\nNumber of sequences downloaded from GenBank :", nrow(lc))
    # cat("\nNumber of sequences selected                :", 1)
    # cat("\nPresent in final alignment                  :", present, "\n\n")
    cat("\n")
    print(head(lc))
    cat("\n")
  }
  
  
  dbDisconnect(conn)  # must be outside of FOR-loop!
  invisible(acc)
}