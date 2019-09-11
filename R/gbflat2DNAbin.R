## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-12-18)

#' @title Read GenBank Flatfiles
#' @description Read a GenBank flatfile, filter it contents according to
#'   organism name and sequence length and return as an object of class
#'   \code{\link{DNAbin}.}
#' @param x A vector of mode \code{"character"} containg the contents of the
#'   GenBank flatfile.
#' @param taxa A vector of mode \code{"character"} giving the set of taxon names
#'   that will be selected from the availbale taxon names (ORGANISM entry in
#'   flatfile).
#' @importFrom ape as.DNAbin
#' @export


gbflat2DNAbin <- function(x, taxa){
  taxa_available <- gsub("(^.+ORGANISM[[:space:]]*)([[:print:]]+)(.+$)", "\\2", x)
  id <- taxa_available %in% taxa
  x <- x[id]
  taxa_available <- taxa_available[id]
  acc <- gsub("(^.+ACCESSION[[:space:]]*)(.+)(\\nVERSION.+$)", "\\2", x)
  dna <- gsub("(^.+ORIGIN[[:space:]]*)(.+)(//$)", "\\2", x)
  dna <- gsub("[[:digit:]]|[[:space:]]|\\n", "", dna)
  
  ## Hier Längenfilter einbauen
  
  dna <- strsplit(dna, "")
  names(dna) <- gsub(" ", "_", paste(taxa_available, acc))
  as.DNAbin(dna)
  
  # names(dna) <- taxa_available
  # write.fas(dna, file = "xxx")
}