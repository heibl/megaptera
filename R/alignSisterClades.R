## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-10-18)

#' @importFrom ips mafft.merge
#' @export

alignSisterClades <- function(tp, seqs, megProj, thread = -1){
   
  seqs <- seqs[as.character(tp)]
  mafft.merge(seqs, exec = megProj@align.exe, 
                      thread = thread)
}