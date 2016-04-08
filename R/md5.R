md5 <- function(x){
  
  write(x, "MD5.txt")
  x <- md5sum("MD5.txt")
  unlink("MD%.txt")
  x
}