## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2016-01-20)

dbUpdateTaxonomy <- function(megapteraProj, taxonomy, tag){
  
  conn2 <- dbconnect(megapteraProj)
  
  ## write data to taxonomy
  ## ----------------------
  if ( !missing(taxonomy) ){
    ## enforce synonym column
    if ( !"synonym" %in% names(taxonomy) ){
      taxonomy <- data.frame(taxonomy, synonym = "-")
    }
    if ( !"tag" %in% names(taxonomy) ){
      if ( missing(tag) ) tag <- "-"
      taxonomy <- data.frame(taxonomy, tag = tag)
    }
    
    taxonomy <- sqlTaxonomyHeader(taxonomy)
    taxonomy$spec <- gsub(" ", "_", taxonomy$spec)
    
    if ( !dbExistsTable(conn2, "taxonomy") ){
      ## create new taxonomy table
      ## -------------------------
      dbWriteTable(conn2, "taxonomy", taxonomy, row.names = FALSE)
      sql <- "ALTER TABLE taxonomy ADD PRIMARY KEY (spec)"
      dbSendQuery(conn2, sql)
      
    } else {
      ## update existing taxonomy table
      ## ------------------------------
      present <- dbReadTable(conn2, "taxonomy")
      
      ## existing table is empty
      ## -----------------------
      if ( nrow(present) == 0 ){

        dbRemoveTable(conn2, "taxonomy")
        dbWriteTable(conn2, "taxonomy", taxonomy, row.names = FALSE)
        sql <- "ALTER TABLE taxonomy ADD PRIMARY KEY (spec)"
        dbSendQuery(conn2, sql)
        
      } else {
        
        ## drop ranks higher than already existing in database
        ## ---------------------------------------------------
        id <- which(!is.na(match(colnames(taxonomy), colnames(present))))[1]
        taxonomy <- taxonomy[, id:ncol(taxonomy)]
        
        ## add columns missing in taxonomy
        ## -------------------------------
        id <- match(colnames(present), colnames(taxonomy))
        taxonomy <- t((t(taxonomy)[id, ]))
        id <- which(is.na(id))
        taxonomy[, id] <- "-"
        colnames(taxonomy)[id] <- colnames(present)[id]
        
        taxonomy <- data.frame(taxonomy)
        dbWriteTable(conn2, "taxonomy", 
                     taxonomy[!taxonomy$spec %in% present$spec, ], 
                     row.names = FALSE, append = TRUE)
      }
    } 
  }
  
  ## check consistency of taxonomy
  
  ## 1. mysterious wrong-genus-error
  ## -------------------------------
  tax <- dbReadTable(conn2, "taxonomy")
  test <- tax$gen != strip.spec(tax$spec)
  if ( any(test) ){
    test <- paste("UPDATE taxonomy",
                  "SET", sql.wrap(strip.spec(tax$spec[test]), term = "gen", BOOL = NULL),
                  "WHERE", sql.wrap(tax$spec[test], term = "spec", BOOL = NULL))
    lapply(test, dbSendQuery, conn = conn2)
  }
  
  ## 2. add species from acc_* tables
  ## --------------------------------
  tabnames <- dbTableNames(megapteraProj, "acc")
  if ( length(tabnames) > 0 ){
    
    ## unset 'excluded (unclassified)'
    SQL <- paste("UPDATE", tabnames,
                 "SET status = 'raw'
                 WHERE status = 'excluded (unclassified)'")
    dbSendQuery(conn2, SQL)
    
    ## tx: taxa in acc_* but not in taxonomy
    ## -------------------------------------
    tx <- paste("SELECT taxon FROM", tabnames,
                "WHERE status !~ 'excluded|too'",
                collapse = " UNION ")
    tx <- paste("(", tx, ") EXCEPT ", 
                "SELECT spec FROM taxonomy",
                sep = "")
    tx <- dbGetQuery(conn2, tx)$taxon
    if ( length(tx) > 0 ){
      gen <- strip.spec(tx)
      s <- paste("SELECT * FROM taxonomy WHERE", 
                 sql.wrap(gen, term = "gen"))
      s <- dbGetQuery(conn2, s)
      if ( nrow(s) == 0 ){ 
        exclude <- tx
      } else {
        id <- gen %in% s$gen
        exclude <- tx[!id]
        
        ## add species to taxonomy table:
        add <- tx[id]
        gen <- gen[id]
        spec.id <- which(colnames(s) == "spec")
        s <- unique(s[, 1:(spec.id - 1)])
        s <- data.frame(s[match(gen, s$gen), ],
                        spec = add, 
                        synonym = "-",
                        tag = "ingroup (NCBI)") ## could also be outgroup!!!
        s <- s[!is.na(s$gen), ]
        dbWriteTable(conn2, "taxonomy", s, 
                     append = TRUE, row.names = FALSE)
        cat("\n", paste(add, collapse = ", "), "added")
      }
      ## mark species as 'excluded (unclassified)'
      ## -----------------------------------------
      if ( length(exclude) > 0 ){
        SQL <- paste("UPDATE", tabnames,
                     "SET status='excluded (unclassified)'")
        SQL <- paste(SQL, "WHERE", sql.wrap(exclude, term = "taxon"))
        lapply(SQL, dbSendQuery, conn = conn2)
        cat("\n", paste(exclude, collapse = ", "), "excluded")
      } 
    }
  }
  
  ## 3. add higher ranks 
  ## -------------------
  hr <- names(tax)[grep("ord", names(tax))[1]]
  fam <- paste("SELECT DISTINCT fam",
               "FROM taxonomy", 
               "WHERE", wrapSQL("-", term = hr, operator = "="), 
               "OR", hr, "IS NULL") 
  fam <- dbGetQuery(conn2, fam)$fam
  
  if ( length(fam) > 0 ){
    cat("\nupdating higher ranks of ")
    for ( i in fam ){
      cat("\n -", i)
      s <- paste("SELECT * FROM taxonomy",
                 "WHERE", wrapSQL(i, term = "fam"))
      s <- dbGetQuery(conn2, s)
      if ( nrow(s) == 1 ){
        cat(": no information")
        next
      }
    
      id <- 1:(which(names(s) == "fam") - 1)
      ## row index where all entries in colums 'id'
      ## are neither NA nor "-
      rid <- is.na(s[, id]) | s[, id] == "-"
      s <- unique(s[!apply(rid, 1, all), id])
      s <- paste(names(s),
                 paste( "'", unlist(s), "'", sep = ""),
                 sep = "=", collapse = ", ")
      s <- paste("UPDATE taxonomy",
                 "SET", s,
                 "WHERE", wrapSQL(i, term = "fam"))
      dbSendQuery(conn2, s)
      cat(": done")
    }
  }
  
  ## 4. delete undeterminet species
  ##    (same set of tokens in stepB + stepBX)
  ## -----------------------------------------
  indet <- c("_sp[.]?([_-]|$)", # Amanita_sp Amanita_sp. Amanita_sp_xxx Amanita_sp._xxx Amanita_sp-53
             "_cf[.]", 
             "_aff[.]", 
             "hybrid(_sp.+)?$", # Juniperus_hybrid Juniperus_hybrid_sp._LO-2009
             "Group$",
             "cultivar$",
             "environmental", # environmental_sample
             "^fungal",
             "uncultured",
             "unknown",
             ".[[:upper:]]",
             "^[[:lower:]]")
  if ( !megapteraProj@taxon@hybrids ){
    indet <- union(indet, "_x_|^x_")
  }
  indet <- paste(indet, collapse = "|")
  indet <- paste("DELETE FROM taxonomy", 
                 "WHERE", wrapSQL(indet, term = "spec", 
                                  boolean = NULL))
  dbSendQuery(conn2, indet)
  # singles quotes are escaped by single quotes in pgSQL!
  dbSendQuery(conn2, "DELETE FROM taxonomy WHERE spec~''''")
  dbDisconnect(conn2)
}

