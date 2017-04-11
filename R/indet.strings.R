## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-01-09)

#' @export

indet.strings <- function(hybrids = TRUE, collapse = FALSE, SQL = FALSE){
  
  obj <- c("( |_)sp[.]?( |-|_|$)", # Amanita_sp  Amanita_sp_xxx
                                   # Amanita_sp. Amanita_sp._xxx Amanita_sp-53
           "spec$",
           "( |_)n[.]sp[.]", # Hydropsyche_n.sp._2006031401
           "( |_)cf[.]", 
           "( |_)aff[.]", 
           # Amylosporus_sp._'succulentus', Limenitis_hybrid_form_'rubidus'
           # singles quotes are escaped by single quotes in pgSQL!
           # e.g. WHERE spec_ncbi~''''
           "'$",
           "hybrid(_sp.+)?$", # Juniperus_hybrid, Juniperus_hybrid_sp._LO-2009, BUT NOT hybridus
           "Group$",
           "cultivar$",
           "environmental", # environmental_sample
           "^fungal",
           "uncultured",
           "unknown",
           ".[[:upper:]]",
           "^[[:lower:]]") 
  
  ## exclude hybrids? (hybrids == FALSE)
  if ( !hybrids ) obj <- c(obj, "_x_", "^x_")
  
  if ( collapse ) obj <- paste(obj, collapse =  "|")
  obj
}