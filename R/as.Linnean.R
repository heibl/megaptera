## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-03-25)

## real names:
## Anigozanthos_humilis_x_Anigozanthos_manglesii
## 'Asterotremella_humicola'
as.Linnean <- function(x){
  
  ## these strings characterize
  ## not fully identified sequences
  ## ------------------------------
  indet <- c("_sp[.]?([_-]|$)", # Amanita_sp Amanita_sp. Amanita_sp_xxx Amanita_sp._xxx Amanita_sp-53
                 "_cf[.]", 
                 "_aff[.]", 
                 "hybrid$", 
                 "Group$",
                 "cultivar$",
                 "environmental$",
                 "^fungal",
                 "uncultured",
                 "unknown",
                 ".[[:upper:]]",
                 "^[[:lower:]]")
  indet <- paste(indet, collapse = "|")
  
  x[grep(indet, x)] <- NA                       
  
  id <- paste("(^[[:upper:]][[:lower:]]+)",  # genus
              "[_| ]",
              "(x?[_| ]?[[:lower:]]+-?[[:lower:]]+)", # epitheton (incl. hybrids)
              #               "([[_| ][var|subsp][.]?]?)", # infraspecific quantifier
              #               "([[_| ][:lower:]+-?[:lower:]+]?)", # subspecies|variety
                            "(.*)", "$", 
              sep = "") 
  x <- gsub(id, "\\1_\\2", x)
  x
}