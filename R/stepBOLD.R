## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-02-20)

#' @export
#' @import DBI

stepBOLD <- function(x, overwrite = TRUE){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  if ( x@locus@kind == "undefined" ) stop("undefined locus not allowed")
  if ( !url.exists("https://eutils.ncbi.nlm.nih.gov") )
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
  
  ## get species names
  ## -----------------
  conn <- dbconnect(x)
  spec <- paste("SELECT taxon",
                "FROM taxonomy",
                "WHERE", wrapSQL("species", "rank", "="),
                "ORDER BY taxon")
  spec <- dbGetQuery(conn, spec)$taxon
  dbDisconnect(conn)
  
  ## sliding window breaks specis names vector in
  ## batches of 100
  ## --------------
  sw <- seq(from = 1, to = length(spec), by = 100)
  sw <- data.frame(from = sw, to = c(sw[-1] -1, length(spec)))
  sw <- apply(sw, 1, function(z, spec) spec[z[1]:z[2]], spec = spec)
  
  ## wrap bold_seq to get DNAbin
  ## ---------------------------
  formatBOLD <- function(spec, marker){
    b <- bold::bold_seq(spec, marker = marker)
    b <- do.call(rbind, b)
    bb <- paste(b[, "name"], b[, "id"])
    bb <- gsub(" ", "_", bb)
    b <- as.list(b[, "sequence"])
    b <- lapply(b, strsplit, split = "")
    b <- lapply(b, unlist)
    b <- lapply(b, tolower)
    names(b) <- bb
    as.DNAbin(b)
  }
  
  ## download and format sequences
  b <- lapply(sw, formatBOLD, marker = marker)
  b <- do.call(c, b)
  
  ## manage duplicates name + id combinations
  ## where do they come from???
  ## --------------------------
  d <- duplicated(names(b))
  if ( any(d) ){
    names(b)[d] <- paste(names(b)[d], 1:length(which(d)), sep = "-")
    # a <- b[names(b) == "Locusta_migratoria_CYTC5284-12"]
    # a <- mafft(a, exec = x@align.exe)
    # rownames(a) <- paste(rownames(a), 1:nrow(a), sep = "_")
    # write.nex(a, "aaa.nex")
  }
  stepBX(x, b, tag = "bold", overwrite = overwrite)
}