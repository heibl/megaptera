## This code is part of the megaptera package
## © C. Heibl 2014 (2016-04-08)

alignGenus <- function(genus, megProj){
  
  ## PARAMETERS
  ## -----------
  gene <- megProj@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- megProj@taxon@tip.rank
  msa.tab <- paste(tip.rank, gsub("^_", "", gene), sep = "_")
  max.bp <- megProj@params@max.bp
  logfile <- paste(gene, "stepG.log", sep = "-")
  align.exe <- megProj@align.exe
  
  ## read sequences from database
  ## ----------------------------
  slog("\n   -", genus, file = logfile)
  if ( tip.rank == "spec" ) genus <- paste(genus, "_", sep = "")
  ## note the factor in max.bp: this is because consensus sequences
  ## can be longer than the longest of its input sequences
  seqs <- dbReadDNA(megProj, tab.name = msa.tab, 
                    taxon = genus, regex = TRUE,
                    max.bp = max.bp * 1.5, 
                    blocks = "ignore")
  
  n <- ifelse(is.list(seqs), length(seqs), nrow(seqs))
  aligned <- ifelse(megProj@update, FALSE, is.matrix(seqs))
  
  if ( !aligned & n > 1 ) {
    slog(" -- aligning", file = logfile)
    seqs <- del.gaps(seqs)
    ref <- dbReadReference(megProj, gene)[1, ] ## dirty subsetting!
    rownames(ref) <- "REF"
    seqs <- c(seqs, as.list(ref))
    seqs <- mafft(seqs, method = "localpair", maxiterate = 1000,
                  path = align.exe, quiet = TRUE)
    seqs <- seqs[rownames(seqs) != "REF", ]
    seqs <- trimEnds(seqs, 1)
    dbWriteMSA(megProj, dna = seqs, status = "genus-aligned")
  }
  else {
    seqs <- as.matrix(seqs)
    if ( tip.rank == "spec" ) slog(" --", n, "species", file = logfile)  
  }
  seqs
}