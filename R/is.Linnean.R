## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-10-24)

#' @export
 
is.Linnean <- function(x, det = TRUE){
  id <- paste0("(^[[:upper:]][[:lower:]]+)",  # genus
              "(_| )",
              "(x?[_| ]?[[:lower:]]+-?[[:lower:]]+)", # epitheton (incl. hybrids)
              #               "([[_| ][var|subsp][.]?]?)", # infraspecific quantifier
              #               "([[_| ][:lower:]+-?[:lower:]+]?)", # subspecies|variety
              "(.*)", "$") 
  id <- grep(id, x)
  seq_along(x) %in% id
}
