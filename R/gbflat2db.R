## This code is part of the megaptera package
## © C. Heibl 2019 (last update 2021-03-18)

#' @title Read GenBank Flatfiles
#' @description Read a GenBank flatfile, optionally filter its contents according to
#'   organism name and sequence length and return as an object of class
#'   \code{\link{DNAbin}.}
#' @param file A vector of mode \code{"character"} the path to the GenBank flatfile.
#' @param taxa A vector of mode \code{"character"} giving the set of taxon names
#'   that will be selected from the availbale taxon names (ORGANISM entry in
#'   flatfile).
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @importFrom ape as.DNAbin
#' @importFrom crayon %+% bold cyan green magenta silver
#' @export


gbflat2db <- function(file, taxa = NULL, megProj){
  
  slog(silver(file) %+% "\n")
  conn <- file(file, open = "r")
  seqs <- readLines(con = conn)
  close(conn)
  if (length(seqs) < 100) {
    stop("debug me!")
  }
  
  ## Identify individual accessions and extract them
  ## ------------------------------------------------
  end <- which(seqs == "//")
  start <- grep("^LOCUS", seqs[1:100])[1]
  start <- c(start, head(end, -1) + 1)
  index <- matrix(c(start, end), ncol = 2)
  seqs <- apply(index, 1, 
                function(z, i) paste0(z[i[1]:i[2]], collapse = '\n'),
                z = seqs)
  remove(start, end, index)
  
  ## Extract study group's sequences
  ## -------------------------------
  if (length(seqs) == 1){
    
    ## This is a hack to avoid the code breaking at the whole genome of Triticum
    ## turgidum subsp. durum (LT934116), which makes up one single flatfile.
    ## gsub does not seem to be able a vector element of size 1GB [2019-10-11]
    ss <- substr(seqs, 1, 1000)
    taxa_available <- gsub("(^.+ORGANISM[[:space:]]*)([[:print:]]+)(.+$)", "\\2", ss)
  } else {
    taxa_available <- gsub("(^.+ORGANISM[[:space:]]*)([[:print:]]+)(.+$)", "\\2", seqs)
  }
  
  if (!is.null(taxa)){
    id <- taxa_available %in% unlist(taxa)
    if (any(id)){
      seqs <- seqs[id]
      taxa_available <- taxa_available[id]
    } else {
      return(NULL)
    }
    remove(id)
  }
  
  ## Create DNAbin list
  ## ------------------
  acc <- gsub("(^.+ACCESSION[[:space:]]*)(.+)(\\nVERSION.+$)", "\\2", seqs)
  acc <- gsub("\\n", "", acc) ## line break in Olea europaea
  seqs <- gsub("(^.+ORIGIN.+ 1 )(.+)(//$)", "\\2", seqs)
  seqs <- gsub("[[:digit:]]|[[:space:]]|\\n", "", seqs)
  ## Hier Längenfilter einbauen
  
  ## Replace synonyms (and surrogates) by the accepted name according to the
  ## user (i.e. the first element of every element of taxa)
  ## ------------------------------------------------------
  chooseAccepted <- function(spec, taxlist){ 
    id <- unlist(lapply(taxlist, function(a, b) b %in% a, b = spec))
    taxlist[[which(id)]][1]
  }
  taxa_accepted <- sapply(taxa_available, chooseAccepted, taxlist = taxa)
  slog(silver(" > found " %+% magenta$bold(length(taxa_available)) %+% " sequences"))
  
  ## Create data frame and write to database
  ## ---------------------------------------
  seqs <- data.frame(acc, taxon = taxa_accepted, 
                     taxon_source = taxa_available, sequence = seqs,
                     stringsAsFactors = FALSE)
  if (nrow(seqs)) {
    conn <- dbconnect(megProj)
    present <- dbGetQuery(conn, "SELECT acc FROM sequence")
    if (nrow(present)){
      id <- seqs$acc %in% present$acc
      if (length(id)){
        ## Only new sequences will be written to database
        ## ----------------------------------------------
        # slog("\n..", length(which(id)), "duplicates removed ..", file = logfile, 
        #      megProj = megProj)
        seqs <- seqs[!id, ]
        
      }
    }
    slog(silver(" > writing " %+% magenta$bold(nrow(seqs)) %+% " sequences to database"))
    dbWriteTable(conn, "sequence", seqs, row.names = FALSE, append = TRUE)
    # slog("\n..", nrow(seqs), "sequences written to", acc.tab, "", 
    # file = logfile, megProj = megProj)  
    dbDisconnect(conn)
  } 
}
