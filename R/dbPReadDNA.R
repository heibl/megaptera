## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-012-23)

dbPReadDNA <- function(conn, tab.name, taxon, regex = FALSE, 
                      max.bp, min.identity, min.coverage,
                      ignore.excluded = TRUE,
                      masked = FALSE){
  
  close.after <- FALSE
  if ( inherits(conn, "megapteraProj") ){
    conn <- dbconnect(conn)
    close.after <- TRUE
  }
  
  if ( missing(taxon) ) taxon <- ".+"
  otaxon <- taxon
  
  dna <- ifelse(masked, "masked", "dna") 
  
  ## esacape metacharacters in taxon names
  ## -------------------------------------
  if ( !regex ){
    taxon <- gsub(" ", "_", taxon)
    taxon <- gsub("[.]$", "[.]", taxon)
    taxon <- gsub("([(]|[+]|[)])", "[\\1]", taxon)
    taxon <- gsub("'", ".", taxon) # e.g.Acorus_sp._'Mt._Emei'
    taxon <- paste("^", taxon, "$", sep = "")
  }
  
  ## field names in <tab.name>
  ## -------------------------
  cols <- paste("SELECT column_name FROM information_schema.columns WHERE", 
                sql.wrap(tab.name, term = "table_name"), 
                "ORDER by ordinal_position")
  cols <- dbGetQuery(conn, cols)$column_name
  
  if ( "taxon" %in% cols ){
    
    ## retrieve sequences from acc.table
    ## ---------------------------------
    SQL <- paste("SELECT taxon, gi, dna FROM ", tab.name, 
                 " WHERE taxon ~ '", taxon, "'", sep = "")
    if ( !missing(max.bp) ) SQL <- paste(SQL, "AND npos <=", max.bp)
    if ( !missing(min.identity) ) SQL <- paste(SQL, "AND identity >=", min.identity)
    if ( !missing(min.coverage) ) SQL <- paste(SQL, "AND coverage >=", min.coverage)
    if ( ignore.excluded ) SQL <- paste(SQL, "AND status !~ 'excluded|too'")
#     print(SQL)
    seqs <- dbGetQuery(conn, SQL)
    if ( nrow(seqs) == 0) stop("no sequences for '", otaxon, "'")
    seq.names <- paste(seqs$taxon, seqs$gi, sep = "_")

  } else {
    
    ## retrieve sequences from spec.table
    ## ----------------------------------
    SQL <- ifelse(ignore.excluded,
                   " AND status !~ 'excluded'", "") 
    SQL <- paste("SELECT spec, ", dna, " FROM ", tab.name, 
                 " WHERE spec ~ '", taxon, "'", SQL, sep = "")
    if ( !missing(max.bp) ) SQL <- paste(SQL, "AND npos <=", max.bp)
    if ( masked ) SQL <- paste(SQL, "AND status = 'masked'") # masking can drop species from alignments!
    seqs <- dbGetQuery(conn, SQL)
    if ( nrow(seqs) == 0 ) stop("no sequences for '", otaxon, "'")
    seq.names <- seqs$spec
  }

  if ( close.after ) dbDisconnect(conn)

  seqs <- as.list(seqs[, dna])
  seqs <- lapply(seqs, s2c)
  names(seqs) <- seq.names
  seqs <- as.DNAbin(seqs)
  
  ## convert to matrix if possible
  ## -----------------------------
  if ( length(unique(sapply(seqs, length))) == 1 ){
    seqs <- as.matrix(seqs)
    #seqs <- deleteEmptyCells(seqs, quiet = TRUE)
  }
  return(seqs)  
}
