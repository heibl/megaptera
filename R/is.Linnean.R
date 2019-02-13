## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-10-17)

#' @export
 
is.Linnean <- function(x, det = TRUE){
  id <- paste0("(^[[:upper:]][[:lower:]-]+)",  # genus (incl. 'Agarico-suber')
              "(_| )",
              "(x?[_| ]?[[:lower:]]+-?[[:lower:]]+)", # epitheton (incl. hybrids)
              #               "([[_| ][var|subsp][.]?]?)", # infraspecific quantifier
              #               "([[_| ][:lower:]+-?[:lower:]+]?)", # subspecies|variety
              "(.*)", "$") 
  id <- grep(id, x)
  seq_along(x) %in% id
}
