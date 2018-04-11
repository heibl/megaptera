## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-10-24)

#' @export
#' @import DBI

stepBX <- function(x, dna, tag = "user-supplied", overwrite = FALSE){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepBX.log")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""), 
       "\nSTEP BX: adding sequences to database\n",
       file = logfile)
  
  ## format sequences
  ## ----------------
  gitax <- splitGiTaxon(names(dna))
  indet <- indet.strings(x@taxon@hybrids, TRUE)
  status <- rep("raw", nrow(gitax))
  indet <- grep(indet, gitax$taxon)
  status[indet] <- "excluded (indet)"
  dna <- as.character(dna)
  dna <- sapply(dna, paste, collapse = "")
  dna <- data.frame(
    gi = gitax$gi,
    taxon = strip.infraspec(gitax$taxon),
    spec_ncbi = gitax$taxon,
    status = status, 
    genom = tag,
    npos = sapply(dna, nchar),
    identity = NA,
    coverage = NA,
    dna = dna,
    stringsAsFactors = FALSE)
    
    
  ## check taxonomy
  
  ## check for duplicates
  ## --------------------
  conn <- dbconnect(x)
  in.db <- paste("SELECT gi || '_' || taxon AS id FROM", acc.tab)
  in.db <- dbGetQuery(conn, in.db)$id
  in.dna <- paste(dna$gi, dna$taxon, sep = "_")
  dup <- in.dna %in% in.db
  if (any(dup)){
    if (overwrite){
      e <- dna[dup, c("gi", "taxon")]
      e <- paste("DELETE FROM", acc.tab,
                 "WHERE", wrapSQL(e$gi, term = "gi", boolean = NULL),
                 "AND", wrapSQL(e$taxon, term = "taxon", boolean = NULL))
      lapply(e, dbSendQuery, conn = conn)
      slog("\n --", length(e), "seqs. will be overwritten",  
           file = logfile)
    } else {
      dna <- dna[!dup, ]
    }
  }
  
  ## write data to pgSQL database
  ## ----------------------------
  if (nrow(dna)) {
    dbWriteTable(conn, acc.tab, dna, 
                 row.names = FALSE, append = TRUE)
    slog("\n --", nrow(dna), "seqs. written to", acc.tab, 
         file = logfile)  
  } else {
    slog("\n -- all sequences already written to", acc.tab, 
         file = logfile)  
  }
  
  ## declare excluded taxa: has been removed upstream to XML2acc (2016-11-03)
  
  ## select sequences if there are > max.gi.per.spec
  dbMaxGIPerSpec(x)
  
  ## create and update relation <taxonomy>
  dbUpdateTaxonomy(x, logfile = logfile) # handle species found in stepB X
  # that are not included in taxonomy table
  
  # summary
  # -------
  total.acc <- dbGetQuery(conn, paste("SELECT count(taxon) FROM", acc.tab))
  det.acc <- dbGetQuery(conn, paste("SELECT count(taxon) FROM", acc.tab, "WHERE status !~'excluded|too'"))
  spec <- dbGetQuery(conn, paste("SELECT taxon, count(taxon) FROM", acc.tab, "WHERE status !~'excluded|too' GROUP BY (taxon)"))
  dbDisconnect(conn)
  one <- length(which(spec$count == 1))
  multiple <- length(which(spec$count > 1))
  det.spec <- one + multiple
  slog(paste("\nNumber of all accessions                              :", total.acc),
       paste("\nNumber of determined accessions                       :", det.acc),
       paste("\nNumber of determined species                          :", det.spec),
       paste("\nNumber of determined species with one accession       : ", one, 
             " (", round(one * 100 / det.spec, 1), "%)", sep = ""),
       paste("\nNumber of determined species with multiple accessions : ", multiple, 
             " (", round(multiple * 100 / det.spec, 1), "%)", sep = ""),
       file = logfile)
  
  
  slog("\n\nSTEP BX finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
}