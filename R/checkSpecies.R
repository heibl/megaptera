## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-05-28)

#' @title Get Information about Species
#' @description Returns information about the 'fate' of a single species along the pipeline.
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @param spec A character string giving the name of a species.
#' @return A data frame, as a side effect a summary message is printed on the screen.
#' @seealso \code{\link{checkMissingSpec}} to create a list of species that were missed/lost during consecutive steps of the pipeline.
#' @import DBI
#' @export

checkSpecies <- function(megProj, spec){
  
  
  conn <- dbconnect(megProj)
  
  ## A: TAXONOMY
  ## -----------
  tax <- dbReadTaxonomy(megProj)
  
  ## if not found, is it a synonym?
  ## ------------------------------
  if (!nrow(tax)){
    tax <- paste("SELECT spec FROM  taxonomy", 
                 "WHERE", wrapSQL(spec, "synonym"))
    tax <- dbGetQuery(conn, tax)$spec
    
    if (!length(tax)){
      cat("\n'", spec, "' not present in 'taxonomy' table\n", sep = "")
      dbDisconnect(conn)
      return()
    } else {
      cat("NOTE:", spec, "is coded as synonym of", tax, "\n")
      spec <- tax
      tax <- paste <- paste("SELECT * FROM  taxonomy", 
                            "WHERE", wrapSQL(tax))
      tax <- dbGetQuery(conn, tax)
    }
  }
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
  # spec <- gsub(" ", "_", spec)
  cat("\n\n***** SEQUENCES *****")
  tab <- dbReadLocus(megProj)
  tab <- tab[rownames(tab) == spec, ]
  if (!nrow(tab)){
    cat("\nno sequences of species '", spec, 
        "' in database\n", sep = "")
    dbDisconnect(conn)
    return()
  }
  
  ## which loci
  ## ----------
  loci <- colnames(tab)[tab > 0]
  loci <- unique(gsub("(gb|sel)_", "", loci))
  if (!length(loci)){
    cat("\nno sequences of species '", spec, 
        "' in database\n", sep = "")
  } else {
    cat("\nnumber of loci:", length(loci), "\n")
  }
  
  ACC <- data.frame()
  for (i in seq_along(loci)){
    
    tabb <- tab[, grep(loci[i], colnames(tab)), drop = FALSE]
    SQL <- paste("SELECT gi, status, genom, npos, identity, coverage",
                 "FROM", paste("acc", loci[i], sep = "_"),
                 "WHERE", wrapSQL(spec, term = "taxon"))
    acc <- dbGetQuery(conn, SQL)
    
    if (ncol(tabb) == 1){
      tabb <- data.frame(tabb, "not yet run", stringsAsFactors = FALSE)
      present <- "no"
    } else {
      present <- grep("selected|masked", tabb[, 2])
      present <- ifelse(length(present), "yes", "no")
    }
    
    ## screen output
    ## -------------
    cat("\n*** LOCUS ", i, ": ", loci[i], " ***", sep = "")
    cat("\nnumber of sequences downloaded from GenBank :", tabb[, 1])
    cat("\nnumber of sequences selected for consensus  :", tabb[, 2])
    cat("\npresent in final alignment                  :", present, "\n\n")
    print(acc)
    cat("\n")
    ACC <- rbind(ACC, cbind(loci[i], acc[-9]))
  }
  dbDisconnect(conn)  # must be outside of FOR-loop!
  invisible(ACC)
}