## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-11-13)

#' @title Truncation of Species Names
#' @description Delete infraspecific names and epithets form species names.
#' @param x A character string giving taxon names of species rank and below.
#' @return A character string giving species names.
#' @details \code{strip.infraspec} tries to guess the separating character used;
#'   both \code{" "} and \code{"_"} are possible, but cannot be mixed as, e.g., 
#'   in \code{"Saxicola torquata_axillaris"}.
#' @examples 
#' strip.infraspec(c("Vipera aspis",
#'                   "Vipera aspis aspis",
#'                   "Vipera aspis atra"))
#'      
#' ## separating characters cannot be mixed!
#' spec <- c("Vipera aspis",
#'           "Vipera_aspis_aspis",
#'           "Vipera_aspis atra")                               
#' strip.infraspec(spec) 
#' strip.infraspec(gsub(" ", "_", spec))
#' strip.infraspec(gsub("_", " ", spec))                   
#' @export

strip.infraspec <- function(x){
  
  if (is.factor(x)) x <- levels(x)[x]
  x <- gsub("_x_", "_x-", x) # handle times symbol in hybrids 1
  ## determine separating character
  ## ------------------------------
  empty <- length(grep(" ", x))
  underscore <- length(grep("_", x))
  if (empty & !underscore) sepchar <- " "
  if (!empty & underscore) sepchar <- "_"
  if (empty & underscore) {
    xx <- strsplit(x, " ")
    xx <- sapply(xx, function(z) z[1])
    sepchar <- ifelse(length(grep("_", xx)), "_", " ")
  }
  x <- strsplit(x, sepchar)
  x <- sapply(x, function(x, sc) paste(x[1:2], collapse = sc), sc = sepchar)
  x <- gsub("_x-", "_x_", x) # handle times symbol in hybrids 2
  x
}