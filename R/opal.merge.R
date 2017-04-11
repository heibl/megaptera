## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-02-01)

#' @export

opal.merge <- function(subMSA, merge.exe, mem = "16G"){
  
  ## write imput files, then clear from memory
  fns <- tempfile(pattern = c("in", "in2", "out"), fileext = ".fas")
  write.fas(subMSA[[1]], fns[1])
  write.fas(subMSA[[2]], fns[2])
  remove(subMSA) ## save memory
  
  ## execute opal
  system2(command = merge.exe, 
          args = c(paste("--mem", mem), 
                   paste("--in", fns[1]),
                   paste("--in2", fns[2]),
                   paste("--out", fns[3])
                   )
          )
  
  ## retrieve output, then remove files
  obj <- read.fas(fns[3])
  file.remove(fns)
  obj
}
