## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2017-02-07)

#' @export

dbUpdateTaxonomy <- function(megProj, taxonomy, tag){
  
  conn <- dbconnect(megProj)
  
  ## write data to taxonomy
  ## ----------------------
  if (!missing(taxonomy)){
    ## enforce synonym column
    if ( !"synonym" %in% names(taxonomy) ){
      taxonomy <- data.frame(taxonomy, synonym = "-")
    }
    if (!"tag" %in% names(taxonomy)){
      if (missing(tag)) tag <- "-"
      taxonomy <- data.frame(taxonomy, tag = tag)
    }
    
    taxonomy <- sqlTaxonomyHeader(taxonomy)
    taxonomy$spec <- gsub(" ", "_", taxonomy$spec)
    
    if ( !dbExistsTable(conn, "taxonomy") ){
      ## create new taxonomy table
      ## -------------------------
      dbWriteTable(conn, "taxonomy", taxonomy, row.names = FALSE)
      sql <- "ALTER TABLE taxonomy ADD PRIMARY KEY (spec)"
      dbSendQuery(conn, sql)
      
    } else {
      ## update existing taxonomy table
      ## ------------------------------
      present <- dbReadTable(conn, "taxonomy")
      
      ## existing table is empty
      ## -----------------------
      if ( nrow(present) == 0 ){
        
        dbRemoveTable(conn, "taxonomy")
        dbWriteTable(conn, "taxonomy", taxonomy, row.names = FALSE)
        sql <- "ALTER TABLE taxonomy ADD PRIMARY KEY (spec)"
        dbSendQuery(conn, sql)
        
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
        dbWriteTable(conn, "taxonomy", 
                     taxonomy[!taxonomy$spec %in% present$spec, ], 
                     row.names = FALSE, append = TRUE)
      }
    } 
  }
  
  ## check consistency of taxonomy
  
  
  ## 2. add species from acc_* tables
  ## --------------------------------
  tabnames <- dbTableNames(megProj, "acc")
  if (length(tabnames)){
    
    ## unset 'excluded (unclassified)'
    SQL <- paste("UPDATE", tabnames,
                 "SET status='raw'",
                 "WHERE status='excluded (unclassified)'")
    lapply(SQL, dbSendQuery, conn = conn)
    
    ## tx: taxa in acc_* but not in taxonomy
    ## -------------------------------------
    tx <- paste("SELECT taxon FROM", tabnames,
                "WHERE status !~ 'excluded|too'",
                collapse = " UNION ")
    tx <- paste("(", tx, ") EXCEPT ", 
                "SELECT taxon FROM taxonomy",
                sep = "")
    tx <- dbGetQuery(conn, tx)$taxon
    if (length(tx)){
      gen <- strip.spec(tx)
      s <- paste("SELECT * FROM taxonomy WHERE", 
                 wrapSQL(gen, "taxon", "="))
      s <- dbGetQuery(conn, s)
      if (!nrow(s)){ 
        exclude <- tx
      } else {
        id <- gen %in% s$taxon
        exclude <- tx[!id]
        
        ## add species to taxonomy table:
        add <- tx[id]
      
        max_id <- "SELECT max(parent_id), max(id) FROM taxonomy"
        max_id <- dbGetQuery(conn, max_id)
        max_id <- max(max_id)
        
        s <- data.frame(parent_id = s$id[match(strip.spec(add), s$taxon)],
                        id = max_id + (1:length(add)),
                        taxon = add,
                        rank = "species",
                        stringsAsFactors = FALSE)
        dbWriteTable(conn, "taxonomy", s, 
                     append = TRUE, row.names = FALSE)
        cat("\n", paste(add, collapse = ", "), "added")
      }
      ## mark species as 'excluded (unclassified)'
      ## -----------------------------------------
      if (length(exclude)){
        SQL <- paste("UPDATE", tabnames,
                     "SET status='excluded (unclassified)'")
        SQL <- paste(SQL, "WHERE", sql.wrap(exclude, term = "taxon"))
        lapply(SQL, dbSendQuery, conn = conn)
        cat("\n", paste(exclude, collapse = ", "), "excluded")
      } 
    }
  }
  
  ## 3. add higher ranks 
  ## -------------------
  hr <- names(tax)[1] # hr: highest rank
  fam <- paste("SELECT DISTINCT fam",
               "FROM taxonomy", 
               "WHERE", wrapSQL("-", term = hr, operator = "="), 
               "OR", hr, "IS NULL") 
  fam <- dbGetQuery(conn, fam)$fam
  
  if ( length(fam) > 0 ){
    cat("\nupdating higher ranks of ")
    for ( i in fam ){
      cat("\n -", i)
      s <- paste("SELECT * FROM taxonomy",
                 "WHERE", wrapSQL(i, term = "fam"))
      s <- dbGetQuery(conn, s)
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
      dbSendQuery(conn, s)
      cat(": done")
    }
  }
  
  ## 4. delete undetermined species
  ##    (same set of tokens in stepB + stepBX)
  ## -----------------------------------------
  indet <- indet.strings(x@taxon@hybrids, TRUE, TRUE)
  indet <- paste("DELETE FROM taxonomy", 
                 "WHERE", wrapSQL(indet, term = "spec", 
                                  boolean = NULL))
  dbSendQuery(conn, indet)
  dbDisconnect(conn)
}

