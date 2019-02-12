## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-12-20)

#' @title Extract Genus Name from Linnean Binomial
#' @description Strips the epitheton from a Linnean binomial and return the
#'   genus name.
#' @param x A vector of mode \code{"character"} containing Linnean binomials.
#'   The separating character can be space (\code{" "}) or underscore
#'   (\code{"_"}).
#' @param mode A vector of mode \code{"character"} chosing between two
#'   algorithms: \code{"strsplit"} is fast, but has failed on avery large vector
#'   (> 170000 element); \code{"regex"} is not as fast, but has handled also
#'   large vectors successfully.
#' @return A vector of mode \code{"character"} containing genus names.
#' @seealso \code{\link{strip.infraspec}}
#' @examples
#' strip.spec("Megaptera_novaengeliae")
#'
#' ## separating characters cannot be mixed!
#' spec <- c("Megaptera_novaengeliae", "Megaptera novaeangliae")
#' strip.spec(spec) # does not work
#' strip.spec(gsub("_", " ", spec)) # this works
#' strip.spec(gsub(" ", "_", spec)) # this works
#' @export

strip.spec <- function(x, mode = "strsplit"){
  if (is.factor(x)) x <- levels(x)[x]
  if (mode == "strsplit"){
    sepchar <- ifelse(length(grep("_", x)) != 0, "_", " ")
    x <- strsplit(x, sepchar)
    return(sapply(x, function(x) x[1]))
  } else {
    core <- function(x){
      
      # x <- gsub("[[:blank:]]+", " ", x) ## eliminate double spaces, etc.
      
      ############
      ## Genus
      ############
      
      ## The string will not be accepted as valid genus name
      ## if it does not start with a capital letter (cond_1), 
      ## or if it start with "Unknown" (cond_2) 
      # cond_1 <- !length(grep("^[[:upper:]]", x))
      # cond_2 <- length(grep("^Unknown", x))
      # if (cond_1 | cond_2){ 
      #   return(out)
      # }
      ## Note the hyphen in the search string! Is e.g. for genus "Agarico-suber"
      gsub("(^[[:upper:]][[:lower:]-]+)( |_)(.*$)", "\\1", x)
    }
    sapply(x, core)
  }
}

