## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-06-12)

## when searching/reading x should be interpretable both as literal or regular expression
## when writing, x is literal and the operator must be the equal sign

#' @export

wrapSQL <- function(x, term = "taxon", operator = "~", boolean = "OR", by = 500) {
  
  if (!is.numeric(x)) {
    
    ## singles quotes are escaped by single quotes in pgSQL!
    x <- gsub("'", "''", x) # e.g. "Gigantochloa_sp._'daluoensis'"
    
#     if ( !literal ){
#       ## convert meta characters to literals
#       x <- gsub("([(]|[)]|[+]|[.])", "[\\1]", x) # e.g. "Clivia_sp._RHA+CA_7b", "Oreobolus_sp._1_(Laegaard_70382)
#       
#       ## add start end end meta character
#       if ( !literal ) x <- paste("^", x, "$", sep = "")
#     }
      x <- paste0("'", x, "'") # REGEX only needs quotation
  }
  
  if (!is.null(term)) x <- paste0(term, operator, x)
  if (!is.null(boolean)){
    id <- seq(from = 1, to = length(x), by = 500)
    id <- paste(id, c(id[-1] - 1, length(x)), sep = ":")
    id <- lapply(id, function(obj) eval(parse(text = obj)))
    x <- sapply(id, function(obj, id) paste(obj[id], collapse = paste("", boolean, "")), obj = x)
  }
  x
}
