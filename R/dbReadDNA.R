## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-06-07)

#' @export
#' @import DBI
#' @importFrom seqinr s2c

dbReadDNA <- function(x, tab.name, taxon, regex = TRUE, 
                      max.bp, min.identity, min.coverage,
                      ignore.excluded = TRUE, subtree = FALSE,
                      blocks = "ignore", masked = FALSE){
  
  if (missing(taxon)) taxon <- ".+"
  otaxon <- taxon
  dna <- ifelse(masked, "masked", "dna")
  blocks <- match.arg(blocks, c("ignore", "split", "concatenate"))
  
  ## esacape metacharacters in taxon names
  ## -------------------------------------
  if (!regex){
    # taxon <- gsub(" ", "_", taxon)
    taxon <- gsub("[.]$", "[.]", taxon)
    taxon <- gsub("([(]|[+]|[)])", "[\\1]", taxon)
    taxon <- gsub("'", ".", taxon) # e.g.Acorus_sp._'Mt._Emei'
    taxon <- paste("^", taxon, "$", sep = "")
  }
  
  ## default: return taxon species.* or genus.* table
  ## ------------------------------------------------
  if (missing(tab.name)) {
    tab.name <- paste(x@taxon@tip.rank, x@locus@sql, sep = "_")
    tab.name <- gsub("__", "_", tab.name)
  }
  if (tab.name == "acc") tab.name <- paste("acc", x@locus@sql, sep = "_")
  if (tab.name == x@taxon@tip.rank) tab.name <- paste(x@taxon@tip.rank, x@locus@sql, sep = "_")
  
  
  ## field names in <tab.name>
  ## -------------------------
  conn <- dbconnect(x)
  cols <- paste("SELECT column_name FROM information_schema.columns WHERE", 
                wrapSQL(tab.name, "table_name", "="), 
                "ORDER BY ordinal_position")
  cols <- dbGetQuery(conn, cols)$column_name
  
  if ("gi" %in% cols){
    
    ## retrieve sequences from acc.table
    ## ---------------------------------
    if (missing(max.bp)) max.bp <- x@params@max.bp
    SQL <- paste("SELECT taxon, gi, dna FROM ", tab.name, 
                 " WHERE taxon ~ '", taxon, "'", sep = "")
    SQL <- paste(SQL, "AND npos <=", max.bp)
    if (!missing(min.identity)) SQL <- paste(SQL, "AND identity >=", min.identity)
    if (!missing(min.coverage)) SQL <- paste(SQL, "AND coverage >=", min.coverage)
    if (ignore.excluded) SQL <- paste(SQL, "AND status !~ 'excluded|too'")
#     print(SQL)
    seqs <- dbGetQuery(conn, SQL)
    if (nrow(seqs) == 0) {
      dbDisconnect(conn)
      warning("no sequences for '", otaxon, "'")
      return(NULL)
    }
    seqnames <- paste(seqs$taxon, seqs$gi)
    seqnames <- gsub(" ", "_", seqnames) # enforce underscore (2017-01-12)

  } else {
    
    ## retrieve sequences from species.table or genus.table
    ## ----------------------------------------------------
    tip.rank <- x@taxon@tip.rank
    tr.sql <- "regexp_replace(taxon, ' ', '_') AS taxon"
    SQL <- ifelse(ignore.excluded,
                   "AND status !~ 'excluded'", "")
    if (!missing(max.bp)) {
      SQL <- paste(paste("AND npos <=", max.bp), SQL)
    }
    if (subtree){
      SQL <- paste("WHERE", wrapSQL(taxon, "subtree"), SQL)
    } else {
      SQL <- paste("WHERE", wrapSQL(taxon, "taxon"), SQL)
    }
    SQL <- paste("SELECT", paste(tr.sql, "status", dna, sep = ", "),
                 "FROM", tab.name, SQL)
    # if ( masked ) SQL <- paste(SQL, "AND status ~ 'masked'") # masking can drop species from alignments!
    seqs <- dbGetQuery(conn, SQL)
    if (!nrow(seqs)) {
      dbDisconnect(conn)
      warning("no sequences for '", otaxon, "'")
      return(NULL)
    }
    na <- is.na(seqs[, 3])
    if (all(na)){
      dbDisconnect(conn)
      return(NULL)
    }
    if (any(na)) {
      nodata <- sort(seqs$taxon[na])
      l <- length(nodata)
      if ( l > 12 ){
        nodata <- c(head(nodata), paste("... [", l, "species in total]"))
      }
      warning(paste(nodata, collapse = ", "), " removed")
      seqs <- seqs[!na, ]
    }
    seqnames <- seqs[, "taxon"]
  }
  dbDisconnect(conn)

  ## convert dataframe to list of class DNAbin
  block <- seqs[, 1:2]
  seqs <- as.list(seqs[, dna])
  seqs <- lapply(seqs, seqinr::s2c)
  names(seqs) <- seqnames
  seqs <- as.DNAbin(seqs)
  
  ## split into blocks (if necessary)
  ## --------------------------------
  if (length(grep("block 2", block)) & blocks %in% c("split", "concatenate")){
    
    block <- split(block$taxon, f = block$status)
    seqs <- lapply(block, 
                   function(a, tips) deleteEmptyCells(as.matrix(a[tips]), 
                                                      quiet = TRUE),
                   a = seqs)
    if (blocks == "concatenate"){
      seqs <- do.call(cbind.DNAbin, c(seqs, fill.with.gaps = TRUE))
    }
  } else {
  
    ## if single block: convert to matrix if possible
    ## ----------------------------------------------
    if (length(unique(sapply(seqs, length))) == 1){
      seqs <- as.matrix(seqs)
      seqs <- deleteEmptyCells(seqs, quiet = TRUE)
    }
  }
  return(seqs)  
}
