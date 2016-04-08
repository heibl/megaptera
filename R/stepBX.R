## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-01-20)

stepBX <- function(x, dna, tag = "user-supplied", overwrite = FALSE){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  
  ## DEFINITIONS
  ## -----------
  acc.tab <- paste("acc", x@locus@sql, sep = "_")
  
  ## iniate logfile
  ## --------------
  logfile <- paste(acc.tab, "stepBX.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""), 
       "\nSTEP BX: adding sequences to database\n",
       file = logfile)
  
  ## format sequences
  ## ----------------
  gitax <- splitGiTaxon(names(dna))
  dna <- as.character(dna)
  dna <- sapply(dna, paste, collapse = "")
  dna <- data.frame(
    gi = gitax$gi,
    taxon = gitax$taxon,
    spec_ncbi = gitax$taxon,
    status = "raw", 
    genom = tag,
    npos = sapply(dna, nchar),
    identity = NA,
    coverage = NA,
    dna = dna,
    stringsAsFactors = FALSE)
  
  ## check for duplicates
  ## --------------------
  conn <- dbconnect(x)
  in.db <- paste("SELECT gi || '_' || taxon AS id FROM", acc.tab)
  in.db <- dbGetQuery(conn, in.db)$id
  in.dna <- paste(dna$gi, dna$taxon, sep = "_")
  dup <- in.dna %in% in.db
  if ( any(dup) ){
    if ( overwrite ){
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
  if ( nrow(dna) > 0 ) {
    dbWriteTable(conn, acc.tab, dna, 
                 row.names = FALSE, append = TRUE)
    slog("\n --", nrow(dna), "seqs. written to", acc.tab, 
         file = logfile)  
  } else {
    slog("\n -- all sequences already written to", acc.tab, 
         file = logfile)  
  }
  
  ## declare excluded taxa 
  ## (same set of tokens in stepB + dbUpdateTaxonomy)
  ## ------------------------------------------------
  indet <- c("_sp[.]?([_-]|$)", # Amanita_sp Amanita_sp. Amanita_sp_xxx Amanita_sp._xxx Amanita_sp-53
             "spec$",
             "_cf[.]", 
             "_aff[.]", 
             "hybrid(_sp.+)?$", # Juniperus_hybrid Juniperus_hybrid_sp._LO-2009
             "Group$",
             "cultivar$",
             "environmental", # environmental_sample
             "^fungal",
             "uncultured",
             "unknown",
             ".[[:upper:]]",
             "^[[:lower:]]") 
  if ( !x@taxon@hybrids ){
    indet <- union(indet, "_x_|^x_")
  }
  indet <- paste(indet, collapse = "|")
  indet <- paste("UPDATE", acc.tab, 
                 "SET status = 'excluded (indet)'",
                 "WHERE", wrapSQL(indet, term = "spec_ncbi", 
                                  boolean = NULL))
  dbSendQuery(conn, indet)
  # singles quotes are escaped by single quotes in pgSQL!
  dbSendQuery(conn, paste("UPDATE", acc.tab, 
                          "SET status = 'excluded (indet)'",
                          "WHERE spec_ncbi~''''"))
  
  ## select sequences if there are > max.gi.per.spec
  dbMaxGIPerSpec(x)
  
  ## create and update relation <taxonomy>
  dbUpdateTaxonomy(x) # handle species found in stepB 
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