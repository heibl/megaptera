## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-09-15)

term <- function(organism, kingdom, locus) {
  
  if ( !is.character(organism) )
    stop("organism must be of mode 'character'")
  if ( !inherits(locus, "locus") )
    stop("locus must be of class 'locus'")
  
  ## kingdom added 2016-01-22
  organism <- paste0(organism, "[orgn]+AND+", kingdom, "[orgn]")
  ## eventually move to higher level:
  ##aliases <- paste("\"", locus@aliases, "\"", sep = "")
  aliases <- gsub(" ", "+", locus@aliases)
  
  ## add search fields *CURRENTLY DISABLED*
  if ( FALSE ){
    aliases <- NCBI.wrap(aliases, field = locus@search.fields)
    not <- NCBI.wrap(not, field = locus@search.fields)
  }
  
  ## create URL using 'sgene' object
  ## -------------------------------
  URL <- vector(length = length(organism))
  for ( i in seq_along(organism) ){
    url <- paste0(organism[i], "+AND+", aliases)
    if ( !"not" %in% locus@not ){
      not <- paste0("\"", locus@not, "\"")
      url <- paste0(url, "+NOT+", not)
    }
    if ( length(url) > 1 ){
      url <- paste0("(", url, ")")
      url <- paste(url, collapse = "+OR+")
    }
    URL[i] <- url
  } # end of FOR-loop over i
  if ( length(URL) > 1 ){
    URL <- paste0("(", URL, ")")
    URL <- paste(URL, collapse = "+OR+")
  }
  URL
}
