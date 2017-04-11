## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2017-03-28)

#' @export
#' @import DBI

dbWriteMSA <- function(megProj, dna, 
                       status = "raw", subtree = NULL,
                       n, md5, masked = FALSE){
  
  gene <- megProj@locus@sql
  tip.rank <- megProj@taxon@tip.rank
  msa.tab <- paste(tip.rank, gsub("^_", "", gene), sep = "_")
  conn <- dbconnect(megProj)
  
  ## prepare data
  ## ------------
  if (is.matrix(dna)) dna <- as.list(dna)
  dna <- as.character(dna)
  npos <- sapply(dna, length)
  dna <- lapply(dna, seqinr::c2s)
  spec <- gsub("_", " ", names(dna))
  
  if (!masked){
    in.db <- paste("SELECT taxon",
                   "FROM", msa.tab)
    in.db <- dbGetQuery(conn, in.db)
    insert <- rep(TRUE, length(spec))
    if (nrow(in.db)) insert[spec %in% in.db$taxon] <- FALSE
    
    nmd5 <- ifelse(missing(md5), "",
                   paste(wrapSQL(n, "n", "="), ",", 
                         wrapSQL(md5, "md5", "="), ","))
    for (i in seq_along(dna)){
      SQL <- ifelse(insert[i],
                    paste("INSERT INTO", msa.tab, 
                          "(taxon, n, md5, subtree, status, npos, dna, masked)",  
                          "VALUES (", 
                          wrapSQL(spec[i], NULL), ",",
                          wrapSQL(n, NULL), ",",
                          wrapSQL(md5, NULL), ",",
                          wrapSQL(subtree, NULL), ",",
                          wrapSQL(status, NULL), ",",
                          wrapSQL(npos[i], NULL), ",", 
                          wrapSQL(dna[[i]], NULL), 
                          ", NULL)"),
                    paste("UPDATE", msa.tab, 
                          "SET", wrapSQL(subtree, "subtree", "="), ",", 
                          nmd5, 
                          wrapSQL(status, "status", "="), ",", 
                          wrapSQL(npos[i], "npos", "="), ",", 
                          wrapSQL(dna[[i]], "dna", "="), 
                          "WHERE", wrapSQL(spec[i], "taxon")))
      dbSendQuery(conn, SQL)
    }
  } else {
    
    ## masked alignment
    SQL <- paste("UPDATE", msa.tab, 
                 "SET", wrapSQL(dna, "masked", "=", NULL), 
                 "WHERE", wrapSQL(spec, "taxon", "=", NULL))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  dbDisconnect(conn)
}
