## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-07-12)

indet.strings <- function(){
  
  c("_sp[.]?([_-]|$)", # Amanita_sp Amanita_sp. Amanita_sp_xxx Amanita_sp._xxx Amanita_sp-53
    "spec$",
    "_n[.]sp[.]", # Hydropsyche_n.sp._2006031401
    "_cf[.]", 
    "_aff[.]", 
    "hybrid(_sp.+)?$", # Juniperus_hybrid Juniperus_hybrid_sp._LO-2009
    "Group$",
    "cultivar$",
    "environmental", # environmental_sample
    "^fungal",
    "uncultured",
    "unknown",
    ".[[:upper:]]",
    "^[[:lower:]]") 
}