## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-11-21)

#' @export

slog <- function(..., file = "megaptera.log", sep = " ", megProj){
  
  ## debug levels
  ## 0: no output
  ## 1: only screen output
  ## 2: only file output
  ## 3: screen and file
  ## 4: screen and file, save data upon foreseeabel errors
  ## 5: screen and file output, save data upon foreseeable and unforeseeable errors
  ## ------------------------------------------------------------------------------
  if (missing(megProj)){
    debug.level <- 3
  } else {
    debug.level <- megProj@params@debug.level
  }
  
  if (debug.level > 0){
    x <- c(...)
    ## screen output
    if (debug.level %in% c(1, 3, 4, 5)){
      cat(x, file = "", sep = sep)
    }
    ## file output
    if (file != "" & debug.level > 1)
      cat(x, file = file, sep = sep, append = TRUE)
  }
}