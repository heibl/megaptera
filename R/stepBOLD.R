## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-11-02)

## TO DO: 1. enable search of all taxa
##        2. enable search of taxa that have no sequences 
##           or have not been selected

## current status:
## if species list only missing ingroup species will be searched for
## if higher taxon the entire taxon will be searched. [2017-10-18]

#' @importFrom bold bold_tax_name
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
  
  ## check if locus is available and assign correct abbreviation
  ## according to BOLDSYSTEMS (NOTE: works only for my personal order of aliases)
  ## ----------------------------------------------------------------------------
  markerSet <- c("COI-5P", "ITS", "matK", "rbcL")
  names(markerSet) <- c("cox1", "its", "matk", "rbcl")
  id <- names(markerSet) %in% c(x@locus@sql, x@locus@aliases)
  if (!any(id)) {
    ## locus is not available on BOLD
    slog("\n\nSTEP BX finished", file = logfile)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), file = logfile)
    return() 
  }
  marker <- markerSet[id]
  
  ## species present in NCBI taxonomy
  ## --------------------------------
  conn <- dbconnect(x)
  spec <- paste("SELECT taxon",
                "FROM taxonomy",
                "WHERE", wrapSQL("species", "rank", "="),
                "ORDER BY taxon")
  spec <- dbGetQuery(conn, spec)$taxon
  dbDisconnect(conn)
  
  ingroup_queried <- x@taxon@ingroup
  id <- sapply(ingroup_queried, function(a, b)  any(a %in% b), b = spec)
  ingroup_missing <- ingroup_queried[!id]
  ingroup_missing <- sapply(ingroup_missing, head, n = 1)
  slog(length(ingroup_missing), "ingroup species are missing on NCBI",
       file = logfile)
  
  ## check whether these taxa are present on BOLD SYSTEMS
  ## ----------------------------------------------------
  onBOLD <- bold_tax_name(ingroup_missing)
  id <- is.na(onBOLD$taxid)
  if (any(id)){
    stop("implement this warning!")
    ingroup_missing <- ingroup_missing[-id]
  }
  
  ## sliding window breaks specis names vector in
  ## batches of 100
  ## --------------
  sw <- seq(from = 1, to = length(ingroup_missing), by = 100)
  sw <- data.frame(from = sw, to = c(sw[-1] -1, length(ingroup_missing)))
  sw <- apply(sw, 1, function(z, spec) spec[z[1]:z[2]], spec = ingroup_missing)
  if (is.null(dim(sw))) dim(sw) <- c(1, 1)
  
  ## wrap bold_seq to get DNAbin
  ## ---------------------------
  formatBOLD <- function(x, species, marker){
    b <- BOLD2megaptera(species, marker) 
    if (!length(b)) return(NA)
    bb <- paste(b[, "taxon"], b[, "id"])
    bb <- gsub(" ", "_", bb)
    b <- as.list(b[, "sequence"])
    b <- lapply(b, strsplit, split = "")
    b <- lapply(b, unlist)
    b <- lapply(b, tolower)
    names(b) <- bb
    as.DNAbin(b)
  }
  
  ## download and format sequences
  b <- apply(sw, 1, formatBOLD, x = x, marker = marker)
  b <- b[!is.na(b)]
  if (!length(b)){
    slog("\n.. no sequences on BOLD either", file = logfile)
    # dbProgress(x, "step_bold", "success") # not yet possible
    return()
  }
  slog("Found", length(b), "sequences", file = logfile)
  
  ## manage duplicates name + id combinations
  ## where do they come from???
  ## --------------------------
  b <- do.call(c, b)
  d <- duplicated(names(b))
  if (any(d)){
    slog("Removed", length(which(d)), "duplicates", file = logfile)
    names(b)[d] <- paste(names(b)[d], 1:length(which(d)), sep = "-")
    # a <- b[names(b) == "Locusta_migratoria_CYTC5284-12"]
    # a <- mafft(a, exec = x@align.exe)
    # rownames(a) <- paste(rownames(a), 1:nrow(a), sep = "_")
    # write.nex(a, "aaa.nex")
  }
  
  ## remove sequences that are not determined on species level
  ## ---------------------------------------------------------
  taxa_bold <- splitGiTaxon(names(b))$taxon
  id <- is.Linnean(taxa_bold)
  b <- b[id]
  
  ## write into pgSQL database
  ## -------------------------
  stepBX(x, b, tag = "bold", overwrite = overwrite)
}