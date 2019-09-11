## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-04-15)

## TO DO: 1. enable search of all taxa
##        2. enable search of taxa that have no sequences 
##           or have not been selected

## current status:
## if species list only missing ingroup species will be searched for
## if higher taxon the entire taxon will be searched. [2017-10-18]

#' @importFrom bold bold_stats bold_tax_name
#' @import DBI
#' @export

stepBOLD <- function(x, overwrite = TRUE){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("undefined locus not allowed")
  if (!url.exists("https://eutils.ncbi.nlm.nih.gov"))
    stop("internet connection required for stepBOLD")
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", x@locus@sql, "-stepBOLD.log")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", 
             packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""), 
       "\nSTEP BOLD: searching and downloading sequences from BOLDSYSTEMS\n",
       file = logfile)
  
  ## Check if locus is available and assign correct abbreviation
  ## according to BOLDSYSTEMS (NOTE: works only for my personal 
  ## order of aliases)
  ## ------------------------------------------------------------
  markerSet <- c("COI-5P", "ITS", "matK", "rbcL")
  names(markerSet) <- c("cox1", "its", "matk", "rbcl")
  id <- names(markerSet) %in% c(x@locus@sql, x@locus@aliases)
  if (!any(id)) {
    ## locus is not available on BOLD
    slog("\n\nSTEP BOLD finished", file = logfile)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), file = logfile)
    return() 
  }
  marker <- markerSet[id]
  
  ## What set of species should be queried?
  ## -------------------------------------
  queried <- union(x@taxon@ingroup, x@taxon@outgroup)
  names(queried) <- sapply(queried, head, 1)
  retrieved <- dbReadLocus(x)
  retrieved <- retrieved[grep("gb", names(retrieved))]
  retrieved[retrieved > 1] <- 1
  retrieved <- rowSums(retrieved)
  retrieved <- names(retrieved)[retrieved > 0]
  missing <- setdiff(, retrieved)
  
  slog(length(missing), "queried species do not have any sequences",
       file = logfile)
  
  
  ## It seem that BOLD treats synonyms akwardly, eg. Acossus
  ## terebra and Lamellocossus terebra are presented as separate 
  ## species. Therefore all name sare treated equally 
  m <- lapply(missing, function(z) data.frame(z[1], z))
  m <- do.call(rbind, m)
  names(m) <- c("concept", "name")
  
  public_seqs <- function(taxon){
    bold_stats(taxon)$records_with_species_name
  }
  m$public <- sapply(m$name[], public_seqs)
  
  ## 
  not_public <- setdiff(m$concept, m$concept[m$public > 0])
  if (any(id)){
    cat("\nWARNING:", length(not_public), "species without public data:",
        formatSpecList(not_public))
    m <- m[m$public > 0, ]
  }
  
  ## To do: throw error if no sequences are available
  # if (!length(b)){
  #   slog("\n.. no sequences on BOLD either", file = logfile)
  #   # dbProgress(x, "step_bold", "success") # not yet possible
  #   return()
  # }
  
  ## Download and format sequences
  ## -----------------------------
  b <- lapply(m$name, BOLD2megaptera, marker = marker, 
              out.format = "DNAbin")
  b <- do.call(c, b)
  d <- duplicated(names(b))
  if (any(d)){
    slog("Removed", length(which(d)), "duplicates", file = logfile)
    b <- b[!d]
  }
  slog("Found", length(b), "sequences", file = logfile)
  
  ## Replace names with accepted names
  ## ---------------------------------
  bb <- splitGiTaxon(names(b))
  bb$taxon <- strip.infraspec(bb$taxon)
  bb$concept <- as.character(m$concept)[match(bb$taxon, m$name)]
  if (any(is.na(bb$concept))){
    "debug me!"
  }
  bb <- paste(bb$concept, bb$gi)
  bb <- gsub(" ", "_", bb)
  names(b) <- bb
  
  ## Write into pgSQL database
  ## -------------------------
  stepBX(x, b, tag = "bold", overwrite = overwrite)
}