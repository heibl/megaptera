## This code is part of the megaptera package
## © C. Heibl 2014 (2017-03-24)

#' @importFrom ape del.gaps
#' @importFrom ips deleteEmptyCells mafft trimEnds
#' @export

alignGenus <- function(genus, megProj){
  
  ## PARAMETERS
  ## -----------
  gene <- megProj@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- megProj@taxon@tip.rank
  msa.tab <- paste(tip.rank, gsub("^_", "", gene), sep = "_")
  max.bp <- megProj@params@max.bp
  logfile <- paste0("log/", gene, "-stepG.log")
  align.exe <- megProj@align.exe
  
  ## read sequences from database
  ## ----------------------------
  slog("\n   -", genus, file = logfile)
  if (tip.rank == "spec") genus <- paste(genus, "_", sep = "")
  seqs <- dbReadDNA(megProj, tab.name = msa.tab, taxon = genus)
  
  n <- ifelse(is.list(seqs), length(seqs), nrow(seqs))
  aligned <- ifelse(megProj@update, FALSE, is.matrix(seqs))
  
  if (!aligned & n > 1) {
    slog(" -- aligning", file = logfile)
    seqs <- del.gaps(seqs)
    ref <- dbReadReference(megProj, gene)[1, ] ## dirty subsetting!
    rownames(ref) <- "REF"
    seqs <- c(seqs, as.list(ref))
    seqs <- mafft(seqs, method = "localpair", maxiterate = 1000,
                  ep = .123, exec = align.exe, quiet = TRUE)
    seqs <- seqs[rownames(seqs) != "REF", ]
    seqs <- trimEnds(seqs, 1)
    seqs <- deleteEmptyCells(seqs)
    dbWriteMSA(megProj, dna = seqs, status = "genus-aligned")
  } else {
    ## "alignment" of one single species
    seqs <- as.matrix(seqs)
    if (tip.rank == "spec") slog(" --", n, "species", file = logfile)  
  }
  seqs
}