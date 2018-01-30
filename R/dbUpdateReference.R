## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-01-30)

#' @export

dbUpdateReference <- function(conn, gene, ref){
  
  if ( !inherits(ref, "DNAbin") ) 
    stop("object 'ref' is not of class 'DNAbin'")
  if ( !is.matrix(ref) ) 
    stop("object 'ref' is not of matrix")
  ref <- lapply(as.character(as.list(ref)), seqinr::c2s)
  
  check <- paste("SELECT taxon FROM reference WHERE", 
                 sql.wrap(gene, term = "gene"), 
                 "AND (", sql.wrap(names(ref), 
                                 term = "taxon"), ")")
  check <- dbGetQuery(conn = conn, check)
  insert <- names(ref)
  update <- intersect(insert, check$taxon)
  insert <- setdiff(insert, update)
  
  if ( !is.null(insert) & length(insert) > 0 ){
    insert <- ref[insert]
    SQL <- paste("INSERT INTO reference ", 
                 "(gene, taxon, reference) ",
                 "VALUES (",
                 wrapSQL(gene, term = NULL), 
                 ", ",
                 wrapSQL(names(insert), 
                          term = NULL, boolean = NULL), 
                 ", ",
                 wrapSQL(insert, 
                          term = NULL, boolean = NULL), 
                 ")", sep = "")
    lapply(SQL, dbSendQuery, conn = conn)
  }
  if ( !is.null(update) & length(update) > 0 ){
    update <- ref[update]
    SQL <- paste("UPDATE reference", 
                 "SET", sql.wrap(update, term = "reference", BOOL = NULL), #",", sql.wrap(maxdist, term = "maxdist"),
                 "WHERE", 
                 sql.wrap(gene, term = "gene"), 
                 "AND", 
                 sql.wrap(names(update), term = "taxon", BOOL = NULL))
    lapply(SQL, dbSendQuery, conn = conn)
  }
}