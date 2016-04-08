## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-01-09)

slog <- function(..., file = "megaptera.log"){
  
  x <- c(...)
  cat(x, file = "")
  if ( file != "" )
    cat(x, file = file, append = TRUE)
}