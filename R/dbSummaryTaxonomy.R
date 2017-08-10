## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-08-10)

#' @importFrom DBI dbDisconnect
#' @export

dbSummaryTaxonomy <- function(megProj){
  
  conn <- dbconnect(megProj)
  
  ntotal <- dbGetQuery(conn, "SELECT count(taxon) FROM taxonomy")$count
  nspec <- dbGetQuery(conn, "SELECT count(taxon) FROM taxonomy WHERE rank = 'species'")$count
  ngen <- dbGetQuery(conn, "SELECT count(taxon) FROM taxonomy WHERE rank = 'genus'")$count
  nfam <- dbGetQuery(conn, "SELECT count(taxon) FROM taxonomy WHERE rank = 'family'")$count
  nord <- dbGetQuery(conn, "SELECT count(taxon) FROM taxonomy WHERE rank = 'order'")$count
  dbDisconnect(conn)
  
  n <- format(c(nspec, ngen, nfam, nord), justify = "right")
  r <- format(paste("number of",  c("species", "genera", "families", "orders")), justify = "left")
  r <- paste(r, n, sep = ": ")
  r <- paste(r, collapse = "\n")
  
  cat("taxonomy has", ntotal, "entries\n")
  cat(r)
  
}
