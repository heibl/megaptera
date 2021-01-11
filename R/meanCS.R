#' @export

meanCS <- function(a){
  
  iupac <- c(a = 136, c = 40, g = 72, t = 24, 
             r = 192, y = 48, s = 96, w = 144, k = 80, m = 160, 
             b = 112, d = 208, h = 176, v = 224)
  iupac <- as.raw(iupac)
  
  obj <- NULL
  for (i in seq(1:nrow(a))){
    obj <- c(obj, mean(attr(a, "cs")[i, ][a[i,] %in% iupac]))
  }
  obj
}