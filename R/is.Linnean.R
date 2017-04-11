## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-03-25)

#' @export
 
is.Linnean <- function(x, det = TRUE){
  id <- paste("(^[[:upper:]][[:lower:]]+)",  # genus
              "(_| )",
              "(x?[_| ]?[[:lower:]]+-?[[:lower:]]+)", # epitheton (incl. hybrids)
              #               "([[_| ][var|subsp][.]?]?)", # infraspecific quantifier
              #               "([[_| ][:lower:]+-?[:lower:]+]?)", # subspecies|variety
              "(.*)", "$", 
              sep = "") 
  id <- grep(id, x)
  seq_along(x) %in% id
}
