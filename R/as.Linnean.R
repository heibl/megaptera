## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-02-22)

## real names:
## Anigozanthos_humilis_x_Anigozanthos_manglesii
## 'Asterotremella_humicola'

#' @export

as.Linnean <- function(x){
  
  ## these strings characterize
  ## not fully identified sequences
  ## ------------------------------
  indet <- indet.strings(FALSE, TRUE)
  x[grep(indet, x)] <- NA                       
  
  id <- paste("(^[[:upper:]][[:lower:]]+)",  # genus
              "[_| ]",
              "(x?[_| ]?[[:lower:]]+-?[[:lower:]]+)", # epitheton (incl. hybrids)
              #               "([[_| ][var|subsp][.]?]?)", # infraspecific quantifier
              #               "([[_| ][:lower:]+-?[:lower:]+]?)", # subspecies|variety
                            "(.*)", "$", 
              sep = "") 
  x <- gsub(id, "\\1 \\2", x)
  x
}