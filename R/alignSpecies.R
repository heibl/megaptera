## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-03-02)

alignSpecies <- function(megProj, spec){
  
  gene <- megProj@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  logfile <- paste(gene, "stepC.log", sep = "-")
  
  conn <- dbconnect(megProj)
  seqs <- dbReadDNA(megProj, acc.tab, spec, 
                    max.bp = 2 * megProj@params@max.bp, 
                    ignore.excluded = TRUE)
  slog(paste("\n-- ", ifelse(is.list(seqs), length(seqs), nrow(seqs)), 
             " seqs. of ", spec,  sep = ""), file = logfile) 
  seqs <- mafft(seqs, method = "auto", path = megProj@align.exe)
  
  
  ## Vitis vinifera + trnLF:
  ## bad sequences make alignment longer than max.bp
  ## THIS IS VERY DIRTY!
#   if ( ncol(seqs) > megProj@params@max.bp ){
#     d <- dist.dna(seqs, model = "N", pairwise.deletion = TRUE, 
#                   as.matrix = TRUE)
#     n.zero <- apply(d, 2, function(x) length(x[x == 0]))
#     n.zero <- which.max(n.zero)
#     exclude <- rownames(d)[d[, n.zero] > 10]
#     ## for very rugged alignments the above will not work
#     ## so instead try this:
#     if ( length(exclude) == 0 ){
#       exclude <- names(which.max(colMeans(d)))
#     }
#     seqs <- deleteEmptyCells(seqs[!rownames(seqs) %in% exclude, ])
#     
#     exclude <- splitGiTaxon(exclude)[, "gi"]
#     SQL <- paste("UPDATE", acc.tab, 
#                  "SET status = 'excluded (too distant)'",
#                  "WHERE", sql.wrap(exclude, term = "gi"))
#     dbSendQuery(conn, SQL)
#     slog("\n-- NOTE:", length(exclude), "seqs. of", spec, "excluded", file = logfile)
#   }
  dbWriteDNA(conn, acc.tab, seqs, enforce.binomial = FALSE, 
             status = "aligned")
  
  ## clear results from previous runs of stepE
  ## -----------------------------------------
  SQL <- paste("UPDATE", acc.tab, 
               "SET coverage = NULL, identity = NULL",
               "WHERE", wrapSQL(spec, term = "taxon", operator = "="))
  dbSendQuery(conn, SQL)
  dbDisconnect(conn)
}