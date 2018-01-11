## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-11-03)

#' @rdname dbDNA
#' @export

dbWriteDNA <- function(conn, tab.name, dna, 
                       enforce.binomial = TRUE, status){
  
  if ( is.matrix(dna) ) dna <- as.list(dna)
  dna <- as.character(dna)
  len <- sapply(dna, length)
  dna <- lapply(dna, seqinr::c2s)
  
  if ( missing(status) ) status <- "raw"  
  
  txt <- splitGiTaxon(names(dna), enforce.binomial, sep = " ")
  
  ## UPDATE or INSERT?
  ## -----------------
  SQL <- paste("SELECT gi, taxon FROM", tab.name)
  acc <- dbGetQuery(conn, SQL)
  id <- paste(txt[, "gi"], txt[, "taxon"]) %in% paste(acc$gi, acc$taxon)
  
  ## INSERT sequences
  ## ----------------
  if ( any(!id) ){
    SQL <- paste(wrapSQL(txt[!id, 1], term = NULL, boolean = NULL),
                 wrapSQL(txt[!id, 2], term = NULL, boolean = NULL),
                 wrapSQL("user-sequence", term = NULL, boolean = NULL),
                 wrapSQL(status, term = NULL, boolean = NULL),
                 len[!id],
                 wrapSQL(dna[!id], term = NULL, boolean = NULL),
                 sep = ", ")
    SQL <- paste("INSERT INTO", tab.name, 
                 "(gi, taxon, spec_ncbi, status, npos, dna)",  
                 "VALUES (", SQL, ")")
    lapply(SQL, dbSendQuery, conn = conn)
  }
  
  ## UPDATE sequences
  ## ----------------
  if ( any(id) ){
    SQL <- paste(
      "UPDATE", tab.name, 
      "SET", wrapSQL(status, 
                     term = "status", 
                     operator = "="), 
      ",",
      wrapSQL(dna[id], 
              term = "dna", 
              operator = "=", 
              boolean = NULL), 
      ",",
      wrapSQL(len[id], 
              term = "npos", 
              operator = "=", 
              boolean = NULL), 
      "WHERE", wrapSQL(txt[id, 1], 
                       term = "gi", 
                       operator = "=", 
                       boolean = NULL), 
      "AND", wrapSQL(txt[id, 2], 
                     term = "taxon", 
                     operator = "=", 
                     boolean = NULL)
    )
    lapply(SQL, dbSendQuery, conn = conn)
  }
}
