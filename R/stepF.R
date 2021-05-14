## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2021-03-20)

#' @title Step F: Select Sequences and Assemble FASTA file
#' @description In \code{stepF} FASTA files will be assembled selecting all
#'   sequences that passed the quality evalution.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param update Logical, if \code{FALSE}, step F is redone from scratch, if
#'   \code{TRUE}, updating is done only if the data or parameters have changed.
#'   If left empty, \code{update} is taken from \code{\link{dbPars}}.
#' @param min.cov Numeric between 0 and 100 giving the minimum coverage (see
#'   \code{qcovs} in BLASTN) to be selected into the final alignment.
#' @param src A vector of mode \code{"character"} optionally restricting the
#'   selection to sequences of a particular source.
#' @param user.addition A vector of mode \code{"character"} giving the names of
#'   taxa that should be included in the final alignment even if their BLAST
#'   scores are bad. This is mainly intended for data exploration and should not
#'   be used to create alignments for analysis.
#' @seealso \code{\link{megapteraProj}}; \code{\link{stepE}} for the preceeding
#'   and \code{\link{stepMAFFT}} for the subsequent step.
#' @importFrom DBI dbGetQuery dbRemoveTable dbSendQuery
#' @importFrom ips write.fas
#' @export

stepF <- function(x, update, max.evalue = 11, min.cov = 25, 
                  src = "all", src.dogma = "limitation", 
                  push.og = FALSE, user.addition){
  
  # max.evalue = 0; min.cov = 20; src = "BOLD"; src.dogma = "limitation"
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("please define locus")
  
  ## check if previous step has been run
  ## -----------------------------------
  status <- dbProgress(x)
  # if (status$step_e == "pending") {
  #   stop("the previous step has not been called yet")
  # }
  # if (status$step_e == "error") {
  #   stop("the previous step has terminated with an error")
  # }
  # if (status$step_e == "failure") {
  #   slog("\nNo data from upstream available - quitting", file = "")
  #   dbProgress(x, "step_f", "failure")
  #   return()
  # }
  # if (status$step_e == "success") {
  #   dbProgress(x, "step_f", "error")
  # }
  
  ## PARAMETERS
  ## -----------
  if (missing(update)) update <- x@update
  gene <- x@locus@sql
  tip.rank <-x@taxon@tip.rank
  align.exe <- x@align.exe
  max.bp <- x@params@max.bp
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepF.log")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(silver(bold(paste("megaptera", packageDescription("megaptera")$Version)) %+% "\n"
              %+% as.character(Sys.time()) %+%  "\n"
              %+% bold("STEP F") %+% ": select sequences for alignment\n"))
  
  #was: "\nSTEP F: construct", tip.rank, "consensus sequences\n", 
 
  
  ## open database connection
  conn <- dbconnect(x)
  
  
  
  ## TO DO: nachfolgenden Code sinnvoll einbetten!
  if (!update){
    slog(silver("Update: no -> deleting data from MSA tables\n"), file = logfile)
    SQL <- paste("DELETE FROM sequence_selected", 
                 "WHERE", wrapSQL(gene, "locus"))
    dbSendQuery(conn, SQL)
  } 
  
  taxa_e <- paste("SELECT taxon FROM sequence", 
                   "WHERE", wrapSQL(gene, "locus", "="),
                   "AND evalue <=", max.evalue,
                   "AND sstrand = 'plus'")
  taxa_e <- dbGetQuery(conn, taxa_e)
  taxa_c <- paste("SELECT taxon FROM sequence", 
                   "WHERE", wrapSQL(gene, "locus", "="),
                   "AND qcovs >=", min.cov, ## coverage
                   "AND sstrand = 'plus'")
  taxa_c <- dbGetQuery(conn, taxa_c)
  taxa_ec <- paste("SELECT taxon FROM sequence", 
                "WHERE", wrapSQL(gene, "locus", "="),
                "AND evalue <=", max.evalue,
                "AND qcovs >=", min.cov, ## coverage
                "AND sstrand = 'plus'")
  taxa_ec <- dbGetQuery(conn, taxa_ec)
  slog(silver("Maximum E-value: ") %+% magenta$bold(max.evalue) %+% "\n")
  slog(silver("Minimum coverage (qcovs): ") %+% magenta$bold(min.cov) %+% "\n")
  slog(silver(" > E-value alone would pick ") %+% magenta(bold(nrow(taxa_e)) %+% " accessions/"
                                            %+% bold(length(unique(taxa_e$taxon))) %+% " taxa\n"))
  slog(silver(" > coverage alone would pick ") %+% magenta(bold(nrow(taxa_c)) %+% " accessions/"
                                            %+% bold(length(unique(taxa_c$taxon))) %+% " taxa\n"))
  slog(silver(" > E-value and coverage pick ") %+% magenta(bold(nrow(taxa_ec)) %+% " accessions/"
                                            %+% bold(length(unique(taxa_ec$taxon))) %+% " taxa\n"))
  
  ## Select accessions that have an Expect value of at least 'max.evalue'
  ## and a coverage of at least 'min.cov'
  ## Note: sstrand = minus: cannot be rc'd, they are just bad
  ## -----------------------------------------------------------
  slog(silver("Select candidate sequences ... "))
  taxa <- paste("SELECT taxon, acc, length, mismatch, evalue, qcovs, pident, source",
               "FROM sequence", 
               "WHERE", wrapSQL(gene, "locus", "="),
               "AND evalue <=", max.evalue,
               "AND qcovs >=", min.cov, ## coverage
               "AND sstrand = 'plus'",
               "ORDER BY (qcovs)")
  taxa <- dbGetQuery(conn, taxa)
  if (!nrow(taxa)) {
    dbDisconnect(conn)
    # slog(paste("\n.. WARNING: no sequences comply with min.identity=", 
    #            min.identity, " and min.coverage=", min.coverage, sep = ""), 
    #      file = logfile)
    slog("\n.. WARNING: no sequences available", file = logfile)
    dbProgress(x, "step_f", "failure")
    return()
  }
  if (tip.rank == "genus"){
    taxa$taxon <- strip.spec(taxa$taxon)
  }
  slog(green("OK\n"))
  
  ## Restrict or prefer source of sequences
  if (!is.null(source)){
    nant <- function(z){magenta(bold(length(z)) %+% " accessions/" %+% bold(length(unique(z))) %+% " taxa")}
    src_set <- unique(taxa$source)
    src <- match.arg(src, src_set)
    src.dogma <- match.arg(src.dogma, c("limitation", "preference"))
    if (src.dogma == "limitation"){
      slog(silver("Select only sequences from " %+% magenta$bold(src), "\n"))
      taxa1 <- taxa[taxa$source == src, ]
      slog(silver(" > i.e. available are " %+% nant(taxa1$taxon) %+% " instead of " %+% nant(taxa$taxon) %+% "\n"))
    } else {
      ## src.dogma == "preference"
      slog(silver("Prefer", "Select only sequences from " 
                  %+% magenta$bold(src), "\n"))
      taxa1 <- taxa[taxa$source == src, ]
      taxa2 <- taxa[taxa$source != src, ]
      taxa2 <- taxa2[!taxa2$taxon %in% taxa1$taxon, ]
      slog(silver(" > " %+% nant(taxa1$taxon) %+% " come from " %+% bold(src) %+% "\n"))
      slog(silver(" > " %+% nant(taxa2$taxon) %+% " come from " 
                  %+% bold(paste(src_set[src_set != src], collapse = "/")) %+% "\n"))
      taxa1 <- rbind(taxa1, taxa2)
    }
    
    
    ## Make sure outgroup is included no whether the source is
    ## -------------------------------------------------------
    if (push.og){
      og <- outgroup(x, sep = " ")
      if (!length(intersect(taxa1$taxon, og))){
        slog(silver("Try to push in outgroup sequences ... "))
        taxa_og <- taxa[taxa$taxon %in% og, ]
        if (!nrow(taxa_og)) {
          slog(red("failed\n > no outgroup accessions comply with current setting of "  
                   %+% bold(paste("min.cov", min.cov, sep = " = ")) %+%  " and "  
                   %+% bold(paste("max.evalue", max.evalue, sep = " = ")) %+% "\n"))
        } else {
          slog(green("OK\n"))
          slog(silver(" > " %+% nant(taxa_og$taxon) %+% " of the outgroup included by force\n"))
          taxa1 <- rbind(taxa1, taxa_og)
        }
      }
    }
    
    taxa <- taxa1
  }
  
  ## Add Grade sensu Geneious
  ## ------------------------
  grade <- function(qcovs, evalue, pident, type = "nucleotide"){
    type <- match.arg(type, c("nucleotide", "protein"))
    mgident <- ifelse(type == "nucleotide", 50, 25)       ## minGradedIdentity
    
    (.5 * qcovs                                           ## fractionCoverage
      + 25 * max(0, (1 - (evalue / 10^-200)))             ## eValue
      + 25 * max(0, (pident - mgident)/(100 - mgident)))  ## percentIdentity
    
  } 
  taxa$grade <- grade(taxa$qcovs, taxa$evalue, taxa$pident, type = "nuc")
  
  
  chooseAcc <- function(tab, spec){
    tab <- tab[tab$taxon == spec, ]
    # tab <- tab[tab$evalue == min(tab$evalue), ] # lowest E-value and
    # tab <- tab[tab$qcovs == max(tab$qcovs), ] # greatest coverage
    # tab <- tab[tab$length == max(tab$length), ] # longest alignment
    tab <- tab[tab$length == max(tab$length), ]
    tab[1, ]
  }
  chosen_acc <- lapply(unique(taxa$taxon), chooseAcc, tab = taxa)
  chosen_acc <- do.call(rbind, chosen_acc)
  
  
  ## Check for existing sequences and delete them
  ## --------------------------------------------
  in_tab <- paste("SELECT count(taxon)",
                  "FROM sequence_selected", 
                  "WHERE", wrapSQL(gene, "locus", "="))
  in_tab <- dbGetQuery(conn, in_tab)$count
  if (in_tab){
    slog("\nDelete", in_tab, "sequences from previous runs of stepF", 
              file = logfile)
    in_tab <- paste("DELETE FROM sequence_selected", 
                    "WHERE", wrapSQL(gene, "locus", "="))
    dbSendQuery(conn, in_tab)
  }
  
  ## User addition: this is intended to test sequences that would not have
  ## passed automatic selection
  ## --------------------------
  selected <- chosen_acc$acc
  if (!missing(user.addition)){
    selected <- c(selected, user.addition)
  }
  
  
  slog(silver("Select best sequence per taxon ... "))
  SQL <- paste("SELECT acc, taxon, sequence",
               "FROM sequence", 
               "WHERE", wrapSQL(selected, "acc", "=", "OR", by = 500))
  seqs <- lapply(SQL, dbGetQuery, conn = conn)
  seqs <- do.call(rbind, seqs)
  slog(green("OK\n"))
  slog(silver(" > " %+% magenta$bold(nrow(seqs)) %+% " sequences selected\n"), file = logfile)
  r <- range(sapply(seqs$sequence, nchar))
  slog(silver(" > sequence lengths: " %+% magenta$bold(paste0(r[1],"-", r[2])) %+% " bp\n"), file = logfile)
  
  ## Do reverse-complements 
  ## ----------------------
  # if ("minus" %in% chosen_acc$sstrand){
  #   id <- seqs$acc %in% chosen_acc$acc[chosen_acc$sstrand == "minus"]
  #   seqs$sequence[id] <- rcString(seqs$sequence[id])
  # }
  
  ## Write chosen accessions to MSA table
  ## -------------------------------------
  slog(silver("Write sequences to datebase ... "))
  SQL <- paste(wrapSQL(gene, term = NULL, boolean = NULL),
               wrapSQL(seqs$taxon, term = NULL, boolean = NULL),
               wrapSQL(seqs$acc, term = NULL, boolean = NULL),
               "'raw'",
               wrapSQL(seqs$sequence, term = NULL, boolean = NULL),
               sep = ", ") 
  SQL <- paste("INSERT INTO sequence_selected",
               "(locus, taxon, acc, status, sequence)",
               "VALUES (", SQL, ")")
  lapply(SQL, dbSendQuery, conn = conn)
  slog(green("OK\n"))
  
  
  ## write to FASTA files
  ## ---------------------
  # obj <- dbReadSequenceSelected(x, status = "raw")
  obj <- dbReadMSA(x)
  fn <- paste0("msa/", gene, "-", length(obj), ".fas")
  slog(silver("Write sequences to " %+% bold(fn) %+% " ..."))
  write.fas(obj, fn)
  slog(green("OK\n"))
  
  dbDisconnect(conn)
  slog(silver("\nSTEP F finished"), file = logfile)
  td <- Sys.time() - start
  slog(silver(" after " %+% cyan(paste(round(td, 2), attr(td, "units"))) %+% "\n"), file = logfile)
  dbProgress(x, "step_f", "success")
}
