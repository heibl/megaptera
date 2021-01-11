## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-10-30)

#' @title Identify Undetermined Species
#' @description Provides a set of regular expressions (\code{\link{regex}}) that
#'   identify undetermined names as they will appear in the NCBI Taxonomy.
#' @param exclude.hybrids Logical, indicating if hybrids should be also identified by
#'   the set of regular expressions.
#' @param collapse Logical, if \code{TRUE} the regular expressions are collected
#'   in one single characters string using the \code{|} operator.
#' @param SQL Logical, \emph{currently unused}.
#' @export
#' @keywords internal

indet.strings <- function(exclude.hybrids = TRUE, collapse = FALSE, SQL = FALSE){
  
  obj <- c("^[[:upper:]][[:lower:]]+$", # "Luciola"
           "( |_)sp[.]?( |-|_|$)", # Amanita_sp  Amanita_sp_xxx
           # Amanita_sp. Amanita_sp._xxx Amanita_sp-53
           "spec$",
           "( |_)n[.]sp[.]", # Hydropsyche_n.sp._2006031401
           "( |_)nr[.]", # "near", e.g. Onthophagus nr. babirussa
           "( |_)cf[.]", # "confer"
           "( |_)aff[.]", # "affinis",
           "( |_)gen[.]", # "affinis",
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
           # ".[[:upper:]]", matches "Picea_engelmannii_x_Picea_glauca"
           "^[[:lower:]]") 
  
  ## Exclude hybrids? 
  if (exclude.hybrids) obj <- c(obj, "_x_", "^x_")
  
  if (collapse) obj <- paste(obj, collapse =  "|")
  obj
}


## For debugging: Which regular expressions matches?
# id <- sapply(obj, grep, x = "Picea_engelmannii_x_Picea_glauca")
# obj[sapply(id, function(z) length(z) == 1)]
