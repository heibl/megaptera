## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-10-25)

#' @export

## Formatting of species lists (see also megaptera2Rmarkdown)
## ---------------------------
formatSpecList <- function(spec.list, n.element, n.item = 2){
  
  
  frmt <- function(z, n.syn = 3){
    if (length(z) > 1){
      if (length(z) > n.syn + 1){
        z <- z[1:(n.syn + 2)]
        z[n.syn + 2] <- "..."
      } 
      z <- paste(z[1], paste0("(syn. ", paste(z[-1], collapse = ", "), ")"))
    }
    paste("\n-", z)
  }
  
  if (!missing(n.element)){
    l <- length(spec.list)
    if (length(spec.list) > n.element){
      spec.list <- spec.list[1:(n.element + 1)]
      spec.list[n.element + 1] <- paste0("[+ ", l - length(spec.list) , " elements]")
    }
  }
  sapply(spec.list, frmt)
}
  
  
