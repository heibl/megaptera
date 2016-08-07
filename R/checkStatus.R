## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-08-03)

checkStatus <- function(x, locus){
  
  if ( missing(locus) ) locus <- x@locus@sql
  
  conn <- dbconnect(x)
  tabs <- dbTableNames(x)
  locus <- gsub("^([[:digit:]])", "_\\1", locus) ## ex: 12s -> _12s
  acc.tab <- paste("acc", locus, sep = "_")
  acc.tab <- gsub("__", "_", acc.tab)
  msa.tab <- paste( x@taxon@tip.rank, locus, sep = "_")
  msa.tab <- gsub("__", "_", msa.tab)
  
  obj <- rep(FALSE, 9)
  names(obj) <- LETTERS[1:9]
  
  ## stepA
  if ( "taxonomy" %in% tabs ) obj["A"] <- TRUE
  
  ## stepB
  if ( acc.tab %in% tabs ) obj["B"] <- TRUE
  
  ## stepC
  cc <- paste("SELECT taxon AS spec,", 
              "count(taxon) AS n,",
              "min(char_length(dna)) = max(char_length(dna)) AS aligned",
              "FROM", acc.tab, 
              "WHERE status !~ 'excluded'",
              "GROUP BY taxon") 
  cc <- dbGetQuery(conn, cc)
  if ( all(cc$aligned) ) obj["C"] <- TRUE
  
  ## stepD
  dd <- dbReadReference(x, locus)
  if ( is.logical(dd) ){
    dbDisconnect(conn)
    return(obj)
  } else {
    obj["D"] <- TRUE
  }
  
  ## stepE
  ee <- paste("SELECT count(taxon)",
              "FROM", acc.tab, 
              "WHERE status !~ 'excluded|too'",
              "AND npos <=", x@params@max.bp,
              "AND identity IS NULL")
  ee <- dbGetQuery(conn, ee)
  if ( ee$count == 0 ) obj["E"] <- TRUE
  
  ## stepF
  if ( msa.tab %in% tabs ){
    obj["F"] <- TRUE
  } else {
    dbDisconnect(conn)
    return(obj)
  }
  
  ## stepG
  gg <- dbReadDNA(x, msa.tab, 
                  taxon = ".+", regex = TRUE,
                  blocks = "ignore")
  if ( is.matrix(gg) ){
    obj["G"] <- TRUE
  } else {
    dbDisconnect(conn)
    return(obj)
  }
  
  
  ## stepH
  hh <- paste("SELECT count(status)",
        "FROM", msa.tab, 
        "WHERE status ~ 'block'")
  hh <- dbGetQuery(conn, hh)
  if ( hh$count > 0 ) obj["H"] <- TRUE
  
  ## stepI
  ii <- dbReadDNA(x, msa.tab, 
                  taxon = ".+", regex = TRUE,
                  blocks = "ignore",
                  masked = TRUE)
  if ( !is.null(ii) ) obj["I"] <- TRUE
  
  dbDisconnect(conn)
  obj
}