## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-02-22)

checkSpecies <- function(megapteraProj, spec){
  
  spec <- gsub(" ", "_", spec)
  conn <- dbconnect(megapteraProj)
  
  ## taxonomy
  ## --------
  tax <- paste <- paste("SELECT * FROM  taxonomy", 
                        "WHERE", wrapSQL(spec))
  tax <- dbGetQuery(conn, tax)
  
  ## if not found, is it a synonym?
  ## ------------------------------
  if ( nrow(tax) == 0 ){
    tax <- paste("SELECT spec FROM  taxonomy", 
                 "WHERE", wrapSQL(spec, "synonym"))
    tax <- dbGetQuery(conn, tax)$spec
    
    if ( length(tax) == 0 ){
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
  tax <- tax[1, tax[1, ] != "-", drop = TRUE]
  tax <- data.frame(rank = names(tax), name = unlist(tax))
  tax <- format(tax, justify = "right")
  tax <- paste(tax$rank, tax$name, sep = " : ")
  cat("***** TAXONOMY *****", paste("\n", tax, sep = ""))
  
  ## sequences
  ## ---------
  cat("\n\n***** SEQUENCES *****")
  tab <- dbReadLocus(megapteraProj)
  tab <-tab[rownames(tab) == spec, ]
  if ( nrow(tab) == 0 ){
    cat("\nno sequences of species '", spec, 
        "' in database\n", sep = "")
    dbDisconnect(conn)
    return()
  }
  
  ## which loci
  ## ----------
  loci <- colnames(tab)[tab > 0]
  loci <- unique(gsub("(gb|sel)_", "", loci))
  if ( length(loci) == 0 ){
    cat("\nno sequences of species '", spec, 
        "' in database\n", sep = "")
  } else {
    cat("\nnumber of loci:", length(loci), "\n")
  }
  
  ACC <- data.frame()
  for ( i in seq_along(loci) ){
    
    tabb <- tab[, grep(loci[i], colnames(tab))]
    SQL <- paste("SELECT gi, status, genom, npos, identity, coverage",
                 "FROM", paste("acc", loci[i], sep = "_"),
                 "WHERE", wrapSQL(spec, term = "taxon"))
    acc <- dbGetQuery(conn, SQL)
    
    
    present <- grep("selected|masked", tabb[, 2])
    present <- ifelse(length(present) > 0, "yes", "no")
    
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