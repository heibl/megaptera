## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2016-08-11)

dbWriteMSA <- function(megapteraProj, dna, 
                       status = "raw", subtree = NULL,
                       n, md5, masked = FALSE){
  
  gene <- megapteraProj@locus@sql
  tip.rank <- megapteraProj@taxon@tip.rank
  tab.name <- paste(tip.rank, gsub("^_", "", gene), sep = "_")
  conn <- dbconnect(megapteraProj)
  
  ## prepare data
  ## ------------
  if ( is.matrix(dna) ) dna <- as.list(dna)
  dna <- as.character(dna)
  npos <- sapply(dna, length)
  dna <- lapply(dna, c2s)
  spec <- names(dna)
  
  if ( !masked ){
    in.db <- paste("SELECT", tip.rank, 
                   "FROM", tab.name)
    in.db <- dbGetQuery(conn, in.db)
    insert <- rep(TRUE, length(spec))
    if ( nrow(in.db) > 0 ){
      insert[spec %in% in.db[, tip.rank]] <- FALSE
    }
    nmd5 <- ifelse(missing(md5), "",
                   paste(wrapSQL(n, term = "n", operator = "="), ",", 
                         wrapSQL(md5, term = "md5", operator = "="), ","))
    for ( i in seq_along(dna) ){
      SQL <- ifelse(insert[i],
                    paste("INSERT INTO", tab.name, 
                          paste0("(", tip.rank, ", n, md5, subtree, status, npos, dna, masked)"),  
                          "VALUES (", 
                          sql.wrap(spec[i], term = NULL), ",",
                          sql.wrap(n, term = NULL), ",",
                          sql.wrap(subtree, term = NULL), ",",
                          sql.wrap(md5, term = NULL), ",",
                          sql.wrap(status, term = NULL), ",",
                          sql.wrap(npos[i], term = NULL), ",", 
                          sql.wrap(dna[[i]], term = NULL), 
                          ", NULL)"),
                    paste("UPDATE", tab.name, 
                          "SET", 
                          wrapSQL(subtree, term = "subtree", operator = "="), ",", 
                          nmd5, 
                          wrapSQL(status, term = "status", operator = "="), ",", 
                          wrapSQL(npos[i], term = "npos", operator = "="), ",", 
                          wrapSQL(dna[[i]], term = "dna", operator = "="), 
                          "WHERE", wrapSQL(spec[i], term = tip.rank)))
      dbSendQuery(conn, SQL)
    }
  } else {
    
    ## masked alignment
    SQL <- paste("UPDATE", tab.name, 
                 "SET", wrapSQL(dna, term = "masked", operator = "=", boolean = NULL), 
                 "WHERE", wrapSQL(spec, term = tip.rank, operator = "=", boolean = NULL))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  dbDisconnect(conn)
}
