## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-08-11)

opal.merge <- function(subMSA, merge.exe){
  
  ## write imput files, then clear from memory
  fns <- c("in.fas", "in2.fas", "out.fas")
  write.fas(subMSA[[1]], fns[1])
  write.fas(subMSA[[2]], fns[2])
  remove(subMSA) ## save memory
  
  ## execute opal
  system2(command = merge.exe, 
          args = c("--in in.fas",
                   "--in2 in2.fas",
                   "--out out.fas"))

  ## retrieve output, then remove files
  obj <- read.fas("out.fas")
  file.remove(fns)
  obj
}
