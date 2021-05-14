## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2021-03-19)

#' @export

splitGiTaxon <- function(string, white.space = " "){
  
  ## Output: white space or underscore between genus and species?
  ws <- c(" ", "_")
  white.space <- match.arg(white.space, ws) ## white space
  
  binomial_term <- "(^[[:upper:]][[:lower:][:space:][:punct:]]+)"
  sep_term <- "([ _])"
  ID_term <- "([[:upper:][:digit:][:punct:]]+$)"
  reg_exp <- paste0(binomial_term, sep_term, ID_term)
  taxon <- gsub(reg_exp, "\\1", string)
  gi <- gsub(reg_exp, "\\3", string)
  if (all(taxon == gi)){
    return(list(taxon = gsub(ws[!ws %in% white.space],
                             ws[ws %in% white.space],
                             taxon)))
  }
  list(taxon = gsub(ws[!ws %in% white.space],
                    ws[ws %in% white.space],
                    gsub(reg_exp, "\\1", string)), 
       gi = gsub(reg_exp, "\\3", string))
}
