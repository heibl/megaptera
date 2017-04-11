## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-01-09)

#' @export

taxdumpDaughters <- function(conn, taxon, indet = indet.strings()){
  
  # 1 get tax_id for Strigidae
  if (length(grep("[[:alpha:]]", taxon)) == 1){
    pid <- dbGetQuery(conn, paste("SELECT id FROM names WHERE", 
                                 wrapSQL(taxon, "taxon", "=")))$id
  } else {
    pid <- taxon
  }
  
  # 2: get daughter nodes
  did <- dbGetQuery(conn, paste("SELECT id, rank FROM nodes WHERE", 
                                wrapSQL(pid, "parent_id", "=")))
  if (nrow(did) == 0) return(NULL) # no daughters available
  
  # 3: get taxon names of daughter nodes
  SQL <- paste("SELECT id, taxon FROM names", 
               "WHERE (", wrapSQL(as.character(did$id), "id", "="),
               ") AND name_class = 'scientific name'")
  taxa <- lapply(SQL, dbGetQuery, conn = conn)
  taxa <- do.call(rbind, taxa)
  
  obj <- data.frame(pid, taxa, did$rank, TRUE, stringsAsFactors = FALSE)
  names(obj) <- c("pid", "id", "taxon", "rank", "explode")
  
  ## remove *species* that do not correspond to structurally
  ## valid Latin binomials according to 'indet'
  ## ------------------------------------------
  notvalid <- grep(indet.strings(collapse = TRUE), obj$taxon)
  notvalid <- intersect(notvalid, which(obj$rank == "species"))
  if (length(notvalid) > 0){
    message(length(notvalid), " taxa removed")
    obj <- obj[-notvalid, ]
  }
  
  obj
}