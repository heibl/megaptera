## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2017-11-21)

#' @rdname dbTaxonomy
#' @import DBI
#' @export

dbUpdateTaxonomy <- function(megProj, taxonomy){
  
  tip.rank <- megProj@taxon@tip.rank
  conn <- dbconnect(megProj)
  
  ## write data to taxonomy
  ## ----------------------
  if (!missing(taxonomy)){
    
    ## CREATE NEW TAXONOMY TABLE
    ## -------------------------
    if (!dbExistsTable(conn, "taxonomy")){
      
      dbWriteTable(conn, "taxonomy", taxonomy, row.names = FALSE)
      sql <- "ALTER TABLE taxonomy ADD PRIMARY KEY (taxon)"
      dbSendQuery(conn, sql)
      
      ## UPDATE EXISTING TAXONOMY TABLE
      ## ------------------------------
    } else {
      
      ## Remove rows of rank below 'tip.rank'
      ## ------------------------------------
      id <- vector()
      this.id <- taxonomy[taxonomy$rank == tip.rank, "id"]
      gain <- length(this.id)
      while (gain > 0){
        this.id <- taxonomy[taxonomy$parent_id %in% this.id, "id"]
        id <- c(id, this.id)
        gain <- length(this.id)
      }
      taxonomy <- taxonomy[!taxonomy$id %in% id, ]
      
      ## Restrict taxonomy to those tip taxa
      ## that are already present in the database
      ## ----------------------------------------
      present <- dbReadTable(conn, "taxonomy")
      present_tips <- present$taxon[present$rank == tip.rank]
      delete_tips <- setdiff(taxonomy$taxon[taxonomy$rank == tip.rank], 
                             present_tips)
      delete_tips <- taxonomy$id[taxonomy$taxon %in% delete_tips]
      for (i in delete_tips){
        d <- i
        repeat {
          pid <- taxonomy$parent_id[taxonomy$id == d[1]]
          if (length(taxonomy$id[taxonomy$parent_id == pid]) > 1) break ## there is another taxon of same rank!
          d <- c(pid, d)
        }
        taxonomy <- taxonomy[!taxonomy$id %in% d, ]
      }
      
      
      ## Write data frame to taxonomy table of SQL database
      ## --------------------------------------------------
      dbRemoveTable(conn, "taxonomy")
      dbWriteTable(conn, "taxonomy", taxonomy, row.names = FALSE)
      sql <- "ALTER TABLE taxonomy ADD PRIMARY KEY (taxon)"
      dbSendQuery(conn, sql)
    } 
  } else {
  	cat("\nNo input taxonomy table")
  }
  
  ## check consistency of taxonomy
  
  
  ## 2. add species from acc_* tables
  ## --------------------------------
  cat("\nSummarizing sequence names across loci ... ")
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
    tx <- paste0("(", tx, ") EXCEPT ", 
                "SELECT taxon FROM taxonomy")
    tx <- dbGetQuery(conn, tx)$taxon
    if (length(tx)){
    	cat(length(tx), "are missing from taxonomy")
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
        SQL <- paste(SQL, "WHERE", wrapSQL(exclude, "taxon", "="))
        lapply(SQL, dbSendQuery, conn = conn)
       # cat("\n", paste(exclude, collapse = ", "), "excluded")
      } 
    } else {
    	cat("none are missing from taxonomy")
    }
  }
  
  ## 3. add higher ranks 
  ## -------------------
  # hr <- names(taxonomy)[1] # hr: highest rank
  # fam <- paste("SELECT DISTINCT fam",
               # "FROM taxonomy", 
               # "WHERE", wrapSQL("-", term = hr, operator = "="), 
               # "OR", hr, "IS NULL") 
  # fam <- dbGetQuery(conn, fam)$fam
  
  # if (length(fam)){
    # cat("\nupdating higher ranks of ")
    # for (i in fam){
      # cat("\n -", i)
      # s <- paste("SELECT * FROM taxonomy",
                 # "WHERE", wrapSQL(i, term = "fam"))
      # s <- dbGetQuery(conn, s)
      # if (nrow(s) == 1){
        # cat(": no information")
        # next
      # }
      
      # id <- 1:(which(names(s) == "fam") - 1)
      # ## row index where all entries in colums 'id'
      # ## are neither NA nor "-
      # rid <- is.na(s[, id]) | s[, id] == "-"
      # s <- unique(s[!apply(rid, 1, all), id])
      # s <- paste(names(s),
                 # paste( "'", unlist(s), "'", sep = ""),
                 # sep = "=", collapse = ", ")
      # s <- paste("UPDATE taxonomy",
                 # "SET", s,
                 # "WHERE", wrapSQL(i, term = "fam"))
      # dbSendQuery(conn, s)
      # cat(": done")
    # }
  # }
  
  ## 4. delete undetermined species
  ##    (same set of tokens in stepB + stepBX)
  ## -----------------------------------------
  cat("\nChecking for undetermined taxon names ... ")
  indet <- indet.strings(megProj@taxon@hybrids, TRUE, TRUE)
  indet <- paste("SELECT taxon FROM taxonomy", 
                 "WHERE", wrapSQL(indet, "taxon", "~", NULL),
                 "AND", wrapSQL(tip.rank, "rank", "="))
  n_indet <- nrow(dbGetQuery(conn, indet))
  if (n_indet){
  	cat(n_indet, "found and removed")
  	indet <- gsub("SELECT taxon", "DELETE", indet)
  	dbSendQuery(conn, indet)
  } else {
  	cat("none found")
  }
  invisible(dbDisconnect(conn))
}

