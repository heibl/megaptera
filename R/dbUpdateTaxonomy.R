## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2018-02-21)

#' @rdname dbTaxonomy
#' @import DBI
#' @export

dbUpdateTaxonomy <- function(megProj, taxonomy, logfile = ""){
  
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
    slog("\nNo input taxonomy table", file = logfile)
  }
  
  ################################
  ## Check consistency of taxonomy
  ################################
  
  ## 2. add species from acc_* tables
  ## --------------------------------
  slog("\nSummarizing sequence names across loci ... ", file = logfile)
  tabnames <- dbTableNames(megProj, "acc")
  if (length(tabnames)){
    
    ## Unset 'excluded (unclassified)'
    ## -------------------------------
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
                 "SELECT taxon FROM taxonomy ORDER BY taxon")
    tx <- dbGetQuery(conn, tx)$taxon
    slog(length(tx), "are missing from taxonomy", file = logfile)
    
    if (length(tx)){
      
      gen <- strip.spec(tx) ## allow to contain duplicates
      s <- paste("SELECT * FROM taxonomy WHERE", 
                 wrapSQL(unique(gen), "taxon", "="),
                 "AND", wrapSQL("genus", "rank", "="))
      s <- lapply(s, dbGetQuery, conn = conn)
      s <- do.call(rbind, s)
      if (!nrow(s)){ 
        exclude <- tx
      } else {
        id <- gen %in% s$taxon
        exclude <- tx[!id]
        add <- tx[id]
        
        ## Add those species that have at least one
        ## congener in taxonomy (genus is ancherage point)
        ## -----------------------------------------------
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
        slog("\n", paste(add, collapse = ", "), "added", file = logfile)
      }
      
      ## Try to classify species using BOLD taxonomy
      ## -------------------------------------------
      slog("\nQuery BOLD for missing taxonomic information", file = logfile)
      at_bold <- lapply(exclude, boldLineage)
      names(at_bold) <- exclude
      id <- sapply(at_bold, is.null)
      exclude <- exclude[id]
      at_bold <- at_bold[!id]
      success <- sapply(at_bold, taxdumpAddNode, x = megProj)
      exclude <- c(exclude, names(at_bold)[!success])
      
      
      ## The remaining species cannot be inserted without external 
      ## information and will be tagged as 'excluded (unclassified)'
      ## -----------------------------------------------------------
      if (length(exclude)){
        SQL <- paste("UPDATE", tabnames,
                     "SET status='excluded (unclassified)'")
        SQL <- paste(SQL, "WHERE", wrapSQL(exclude, "taxon", "="))
        lapply(SQL, dbSendQuery, conn = conn)
        write(exclude, file = "log/unclassified_species.txt")
        n_exclude <- length(exclude)
        if (n_exclude > 12) exclude <- c(head(exclude), "[...]", tail(exclude))
        slog("\n", n_exclude, " species tagged as 'excluded (unclassified)':",
             paste("\n-", exclude), file = logfile, sep = "")
        slog("\nSee 'log/unclassified_species.txt' for the full list")
      } 
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
  # slog("\nupdating higher ranks of ", file = logfile)
  # for (i in fam){
  # slog("\n -", i, file = logfile)
  # s <- paste("SELECT * FROM taxonomy",
  # "WHERE", wrapSQL(i, term = "fam"))
  # s <- dbGetQuery(conn, s)
  # if (nrow(s) == 1){
  # slog(": no information", file = logfile)
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
  # slog(": done", file = logfile)
  # }
  # }
  
  ## 4. delete undetermined species
  ##    (same set of tokens in stepB + stepBX)
  ## -----------------------------------------
  slog("\nChecking for undetermined taxon names ... ", file = logfile)
  indet <- indet.strings(megProj@taxon@hybrids, TRUE, TRUE)
  indet <- paste("SELECT taxon FROM taxonomy", 
                 "WHERE", wrapSQL(indet, "taxon", "~", NULL),
                 "AND", wrapSQL(tip.rank, "rank", "="))
  n_indet <- nrow(dbGetQuery(conn, indet))
  if (n_indet){
    slog(n_indet, "found and removed", file = logfile)
    indet <- gsub("SELECT taxon", "DELETE", indet)
    dbSendQuery(conn, indet)
  } else {
    slog("none found", file = logfile)
  }
  invisible(dbDisconnect(conn))
}

