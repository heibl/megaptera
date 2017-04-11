## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-02-17)

#' @export
#' @import DBI

checkStatus <- function(x, locus){
  
  if (missing(locus)) locus <- x@locus@sql
  
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
  ## Criterion: a table names taxonomy exists
  ## ----------------------------------------
  if ("taxonomy" %in% tabs) obj["A"] <- TRUE
  
  ## stepB
  ## Criterion: a table of accessions exists
  ## ---------------------------------------
  if (acc.tab %in% tabs) obj["B"] <- TRUE
  
  ## stepC
  ## Criterion: sequences of each species are
  ## of the same length
  ## ------------------
  cc <- paste("SELECT taxon AS spec,", 
              "count(taxon) AS n,",
              "min(char_length(dna)) = max(char_length(dna)) AS aligned",
              "FROM", acc.tab, 
              "WHERE status !~ 'excluded'",
              "GROUP BY taxon") 
  cc <- dbGetQuery(conn, cc)
  if (all(cc$aligned)) obj["C"] <- TRUE
  
  ## stepD
  ## Criterion: reference table exists
  ## ---------------------------------
  dd <- dbReadReference(x, locus)
  if (is.logical(dd)){
    dbDisconnect(conn)
    return(obj)
  } else {
    obj["D"] <- TRUE
  }
  
  ## stepE
  ## Citerion: the identity column in the
  ## accession table is empty
  ## ------------------------
  ee <- paste("SELECT count(taxon)",
              "FROM", acc.tab, 
              "WHERE status !~ 'excluded|too'",
              "AND npos <=", x@params@max.bp,
              "AND identity IS NULL")
  ee <- dbGetQuery(conn, ee)
  if (!ee) obj["E"] <- TRUE
  
  ## stepF
  ## Citerion:
  ## -
  if (msa.tab %in% tabs){
    
    ## check if msa table is empty
    ff <- paste("SELECT count(dna) FROM", msa.tab)
    ff <- dbGetQuery(conn, ff)
    if ( ff$count > 0 ){
      obj["F"] <- TRUE
    } else {
      dbDisconnect(conn)
      return(obj)
    }
    
  } else {
    dbDisconnect(conn)
    return(obj)
  }
  
  ## stepG
  gg <- dbReadDNA(x, msa.tab)
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
  if (hh$count) obj["H"] <- TRUE
  
  ## stepI
  ii <- dbReadDNA(x, msa.tab, masked = TRUE)
  if ( !is.null(ii) ) obj["I"] <- TRUE
  
  dbDisconnect(conn)
  obj
}