## This code is part of the megaptera package
## © C. Heibl 2018 (last update: 2018-01-22)

#' @title Check A Taxonomic Subset
#' @description Checks the availability of a set of species in a megaptera
#'   project.
#' @param x Either an object of class \code{\link{megapteraProj}} or a data
#'   frame as returned by \code{\link{dbReadTaxonomy}}.
#' @param subset A list of species names. The list's elements can list more than
#'   one species; in this case the first name is taken as accepted name and the
#'   others as synonms.
#' @return A data frame with two columns: 
#' \describe{ 
#' \item{\code{taxon}}{the
#'   accepted name of the taxon} 
#' \item{\code{present}}{the status of the taxon,
#'   i.e. present or not}}
#' @seealso \code{\link{checkSpecies}}
#' @export

checkSubset <- function(x, subset){
  
  
  ## Get species names from NCBI taxonomy
  ## ------------------------------------
  ## get taxonomy if necessary (i.e. tax is not a parent-child-table)
  ## --------------------------------------------------------------
  if (inherits(x, "megapteraProj")){
    x <- dbReadTaxonomy(x)
  }
  x <- x$taxon[x$rank == "species"]
  
  ## Do the matching
  ## ---------------
  core <- function(a, b){
    id <- which(a %in% b)
    if (!length(id)) id <- 0
    if (length(id) > 1) id <- -1
    if (id > 2) id <- 2
    id <- as.integer(unique(id))
  }
  id <- sapply(subset, core, b = x)
  res <- vector(mode = "character", length = length(id))
  res[id == 0] <- "no"
  res[id == 1] <- "yes (accepted)"
  res[id == 2] <- "yes (synonym)"
  res[id == -1] <- "yes (accepted + synonym)"
  res <- data.frame(taxon = names(subset),
                    present = res,
                    stringsAsFactors = FALSE)
  
  ## Screen output
  ## -------------
  ntotal <- length(subset)
  z <- c(nmissing = length(id[id == 0]),
         npresent = length(id[id != 0]),
         nacc = length(id[id == 1]),
         nsyn = length(id[id == 2]),
         naccsyn = length(id[id == -1]))
  zz <- round(z/ntotal * 100, digits = 2)
  z <- format(z, justify = "right")
  zz <- format(paste0("(", zz, "%)"), justified = "right")
  zzz <- c("are missing", "are present", "are present as accepted name",
           "are present as synonym", "are present as accepted name and synonym")
  z <- paste(z, zz, zzz)
  cat(paste0(ntotal, " subset taxa matched against ", length(x), " species in database:"),
      paste0("\n", z))
  
  ## Return results
  res
}