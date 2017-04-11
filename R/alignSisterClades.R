## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-11-09)

#' @export

alignSisterClades <- function(tp, seqs, megProj, thread = -1){
   
  seqs <- seqs[as.character(tp)]
  mafft.merge(seqs, mafft.exe = megProj@align.exe, 
                      thread = thread)
}