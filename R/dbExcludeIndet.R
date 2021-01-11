## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-10-30)

dbExcludeIndet <- function(x){
  
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  indet <- c("_sp[.]?([_-]|$)", # Amanita_sp Amanita_sp. Amanita_sp_xxx Amanita_sp._xxx Amanita_sp-53
             "spec$",
             "_cf[.]", 
             "_aff[.]", 
             "hybrid$", 
             "Group$",
             "cultivar$",
             "environmental", # environmental_sample
             "^fungal",
             "uncultured",
             "unknown",
             ".[[:upper:]]",
             "^[[:lower:]]") 
  if (x@taxon@exclude.hybrids){
    indet <- union(indet, "_x_|^x_")
  }
  indet <- paste(indet, collapse = "|")
  conn <- dbconnect(x@db)
  indet <- paste("UPDATE", acc.tab, 
                 "SET status = 'excluded (indet)'",
                 "WHERE", wrapSQL(indet, term = "spec_ncbi", 
                                  boolean = NULL))
  dbSendQuery(conn, indet)
  # singles quotes are escaped by single quotes in pgSQL!
  dbSendQuery(conn, paste("UPDATE", acc.tab, 
                          "SET status = 'excluded (indet)'",
                          "WHERE spec_ncbi~''''"))
  dbDisconnect(conn)
}