## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-08-10)

#' @title Numerical Summary of Taxonomy
#' @description Returns the number of total entries (nodes) in the taxonomy as
#'   well as the number of species, genera, families and orders.
#' @param megProj An object of class
#'   \code{\linkS4class{megapteraProj}}.
#' @seealso \code{\link{ncbiTaxonomy}} for downloading the NCBI taxonomy;
#'   \code{\link{stepA}} for extracting the taxonomic information relevant for a
#'   given megaptera project; \code{\link{dbUpdateTaxonomy}} and
#'   \code{\link{dbReadTaxonomy}} for storing and retrieving taxonomic
#'   information in a megaptera project database.
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
  
  obj <- c(nspec = nspec, ngen = ngen, nfam = nfam, nord = nord)
  n <- format(obj, justify = "right")
  r <- format(paste("number of",  c("species", "genera", "families", "orders")), justify = "left")
  r <- paste(r, n, sep = ": ")
  r <- paste(r, collapse = "\n")
  
  cat("taxonomy has", ntotal, "entries\n")
  cat(r)
  invisible(c(ntotal, obj))
}
