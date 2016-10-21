## This code is part of the megaptera package
## © C. Heibl 2014 (2016-08-11)

alignSubtree <- function(subtree, taxon, megProj){
  
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
  # slog("\n   -", genus, file = logfile)
  # if ( tip.rank == "spec" ) genus <- paste(genus, "_", sep = "")
  taxon <- taxon$taxon[taxon$subtree == subtree]
  taxon <- paste0("^", taxon, "$")
  taxon <- paste(taxon, collapse = "|")
  seqs <- dbReadDNA(megProj, tab.name = msa.tab, 
                    taxon = taxon, regex = TRUE,
                    blocks = "ignore")
  
  n <- ifelse(is.list(seqs), length(seqs), nrow(seqs))
  aligned <- ifelse(megProj@update, FALSE, is.matrix(seqs))
  slog("\n -- aligning subtree", subtree, "of size", n, file = logfile)
  
  if ( !aligned & n > 1 ) {
    
    # seqs <- del.gaps(seqs)
    # ref <- dbReadReference(megProj, gene)[1, ] ## dirty subsetting!
    # rownames(ref) <- "REF"
    # seqs <- c(seqs, as.list(ref))
    seqs <- mafft(seqs, method = "localpair", maxiterate = 1000,
                  path = align.exe, quiet = TRUE)
    # seqs <- seqs[rownames(seqs) != "REF", ]
    seqs <- trimEnds(seqs, 1)
    seqs <- deleteEmptyCells(seqs, quiet = TRUE)
    dbWriteMSA(megProj, dna = seqs, subtree = subtree, 
               status = "subtree-aligned")
  } else {
    ## "alignment" of one single species
    seqs <- as.matrix(seqs)
    dbWriteMSA(megProj, dna = seqs, subtree = subtree, 
               status = "subtree-aligned")
  }
  # seqs
}