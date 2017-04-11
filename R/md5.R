## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-03-28)

#' @export

md5 <- function(x){
  
  fn <- tempfile("md5")
  write(x, fn)
  x <- tools::md5sum(fn)
  unlink(fn)
  x
}