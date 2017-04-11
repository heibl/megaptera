## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-11-17)

#' @export
#' @import DBI

stepPASTAout <- function(x){	
  
  start <- Sys.time()
  quiet = FALSE
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  if ( x@locus@kind == "undefined" ) 
    stop("undefined locus not allowed")
  STATUS <- checkStatus(x)
  if ( !all(STATUS[1:6]) ){
    stop("step", names(STATUS)[min(which(!STATUS))] ,
         " has not been run")
  }
  
  ## PARAMETERS
  ## ----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- x@taxon@tip.rank
  msa.tab <- paste(tip.rank, gsub("^_", "", gene), sep = "_")
  min.n.seq <- x@params@min.n.seq
  align.exe <- x@align.exe
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "stepPASTAout.log")
  if (!quiet & file.exists(logfile)) unlink(logfile)
  if (!quiet)  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
                      paste("\n", Sys.time(), sep = ""),
                      "\nSTEP G: alignment\n", 
                      paste("\n.. locus:", x@locus@sql), file = logfile)
  
  seqs <- paste0(gene, "[[:digit:]][.]marker001[.]", gene, ".aln")
  seqs <- list.files(pattern = seqs)
  seqs <- read.fas(tail(seqs, 1))
  
  ## prune ends of alignment from 'thin tails'
  ## -----------------------------------------
  # seqs <- trimEnds(seqs, min(nrow(seqs), min.n.seq))
  seqs <- deleteGaps(seqs, min.gap = nrow(seqs) - 1)
  seqs <- deleteEmptyCells(seqs, quiet = TRUE)
  
  ## sort alignment taxonomically
  ## ----------------------------
  if (nrow(seqs) > 1){
    gt <- dbReadTaxonomy(x, subset = seqs)
    gt <- tax2tree(gt, tip.rank = tip.rank)
    gt <- ladderize(gt)
    gt <- drop.tip(gt, setdiff(gt$tip.label, rownames(seqs)))
    seqs <- seqs[match(gt$tip.label, rownames(seqs)), ]
  }
  
  ## write to database -- serially or in parallel
  ## --------------------------------------------
  ## open database connection
  conn <- dbconnect(x)
  slog("\n.. write alignment to database ..", file = logfile)
  # if ( cluster.open ){
  #   if ( nrow(seqs) > 1000 ){
  #     id <- seq(from = 1, to = nrow(seqs), by = ceiling(nrow(seqs)/x@params@cpus))
  #     id <- data.frame(from = id, to = c(id[-1] - 1, nrow(seqs)))
  #     aa <- apply(id, 1, function(i, a) a[i[1]:i[2], ], a = seqs)
  #     sfLibrary("seqinr", character.only = TRUE)
  #     sfExport("aa")
  #     sfLapply(x = aa, fun = dbWriteMSA, megapteraProj = x, 
  #              status = "aligned")
  #   } else {
  #     dbWriteMSA(x, dna = seqs, status = "aligned")
  #   }
  #   sfStop()
  # } else {
    dbWriteMSA(x, dna = seqs, status = "aligned")
  # }
  ## This is just a hack to allow for negative matching of regexs
  ## In the future stepG should insert the highest non-saturated taxonomic rank
  dbSendQuery(conn, paste("UPDATE", msa.tab, 
                          "SET status = 'included'", 
                          "WHERE status is NULL"))
  dbDisconnect(conn)
  
  ## write files
  ## -----------
  write.phy(seqs, paste(gene, "phy", sep = "."))
  rownames(seqs) <- gsub("-", "_", rownames(seqs))
  write.nex(seqs, paste(gene, "nex", sep = "."))
  
  ## calculate mean absolute deviation (MAD)
  ## see Smith, Beaulieau, Donoghue (2009)
  this.mad <- round(MAD(seqs), 5)
  slog("\n.. mean absolute deviation:", 
       this.mad, "..", file = logfile)
  
  ## summary
  ## -------
  if (!quiet) {
    slog(
      paste("\n\n--- final alignment of", gene, "---"),
      paste("\nnumber of sequences     :", nrow(seqs)),
      paste("\nnumber of base pairs    :", ncol(seqs)),
      paste("\nmean absolute deviation :", this.mad),
      "\n\nSTEP G finished",
      file = logfile)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), file = logfile)
  }
  invisible(x)
}
