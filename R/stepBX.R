## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2021-03-19)

#' @export
#' @import DBI
#' @importFrom crayon %+% bold cyan magenta silver

stepBX <- function(x, dna, tag = "user-supplied", overwrite = FALSE){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("undefined locus not allowed")
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- "sequence"
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepBX.log")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(silver(bold(paste("\nmegaptera", packageDescription("megaptera")$Version)),
       paste("\n", Sys.time(), "\n", sep = "") %+%
       bold("STEP BX") %+% ": adding sequences to database\n"),
       file = logfile)
  
  ## Sequences must not be aligned
  ## ----------------------------
  slog(silver("Removing gaps from alignment ... "))
  dna <- del.gaps(dna)
  slog(green("OK\n"))
  
  ## Sequences must not be aligned
  ## ----------------------------
  slog(silver("Extracting study group sequences ... "))
  tax <- dbReadTaxonomy(x)
  gitax <- as.data.frame(splitGiTaxon(names(dna), white.space = " "))
  id <- gitax$taxon %in%  tax$taxon[tax$rank == x@taxon@tip.rank]
  dna <- dna[id]
  gitax <- gitax[id, ]
  slog(silver(green("OK\n") %+% " > " %+% magenta$bold(nrow(gitax)) %+% " sequences (of " %+% 
         magenta$bold(length(unique(gitax$taxon))) %+% " taxa) belong to the study group" %+% "\n"))
  
  ## format sequences
  ## ----------------
  slog(silver("Formatting sequences ... "))
  indet <- indet.strings(x@taxon@exclude.hybrids, TRUE)
  status <- rep("raw", nrow(gitax))
  indet <- grep(indet, gitax$taxon)
  status[indet] <- "excluded (indet)"
  dna <- as.character(dna)
  dna <- sapply(dna, paste, collapse = "")
  dna <- data.frame(
    acc = gitax$gi,
    taxon = strip.infraspec(gitax$taxon),
    taxon_source = gitax$taxon,
    locus = gene,
    status = status, 
    sequence = dna,
    stringsAsFactors = FALSE)
  slog(green("OK\n"))
  
  ## Check for duplicates
  ## --------------------
  slog(silver("Checking for duplicates ... "))
  conn <- dbconnect(x)
  in.db <- paste("SELECT acc || '_' || taxon AS id FROM", acc.tab)
  in.db <- dbGetQuery(conn, in.db)$id
  in.dna <- paste(dna$acc, dna$taxon, sep = "_")
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
  slog(green("OK\n"))
  
  ## Write data to pgSQL database
  ## ----------------------------
  slog(silver("Writing sequences to database ... \n"))
  if (nrow(dna)) {
    dbWriteTable(conn, acc.tab, dna, 
                 row.names = FALSE, append = TRUE)
    slog(silver(" > " %+% magenta$bold(nrow(dna)) %+% "seqs. written to"), acc.tab, 
         file = logfile)  
  } else {
    slog(silver(" > all sequences already written to " %+% acc.tab), 
         file = logfile)  
  }
  
  ## declare excluded taxa: has been removed upstream to XML2acc (2016-11-03)
  
  ## select sequences if there are > max.gi.per.spec
  # dbMaxGIPerSpec(x)
  
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
  slog(silver("Number of all accessions                              : "
              %+% magenta$bold(total.acc) %+% "\n"
              %+% "Number of determined accessions                       : "
              %+% magenta$bold(det.acc) %+% "\n"
              %+% "Number of determined species                          : "
              %+% magenta$bold(det.spec) %+% "\n"
              %+% "Number of determined species with one accession       : "
              %+% magenta(bold(one) %+% paste0(" (", round(one * 100 / det.spec, 1), "%)\n"))
              %+% "Number of determined species with multiple accessions : " 
              %+% magenta(bold(multiple) %+% paste0(" (", round(multiple * 100 / det.spec, 1), "%)\n"))),
       file = logfile)
  
  
  slog("\n\nSTEP BX finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
}
