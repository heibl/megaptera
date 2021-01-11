## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-10-30)

#' @title Step F: Select Sequences and Assemble FASTA file
#' @description In \code{stepF} FASTA files will be assembled selecting all
#'   sequences that passed the quality evalution.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param update Logical, if \code{FALSE}, step F is redone from scratch, if
#'   \code{TRUE}, updating is done only if the data or parameters have changed.
#'   If left empty, \code{update} is taken from \code{\link{dbPars}}.
#' @param min.cov Numeric between 0 and 100 giving the minimum coverage (see
#'   \code{qcovs} in BLASTN) to be selected into the final alignment.
#' @param user.addition A vector of mode \code{"character"} giving the names of
#'   taxa that should be included in the final alignment even if their BLAST
#'   scores are bad. This is mainly intended for data exploration and should not
#'   be used to create alignments for analysis.
#' @seealso \code{\link{megapteraProj}}; \code{\link{stepE}} for the preceeding
#'   and \code{\link{stepMAFFT}} for the subsequent step.
#' @importFrom DBI dbGetQuery dbRemoveTable dbSendQuery
#' @importFrom ips write.fas
#' @export

stepF <- function(x, update, min.cov = 25, user.addition){
  
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
  if (file.exists(logfile)) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP F: construct", tip.rank, "consensus sequences\n", 
       file = logfile)
  
  ## open database connection
  conn <- dbconnect(x)
  
  ## TO DO: nachfolgenden Code sinnvoll einbetten!
  if (!update){
    slog("\nUpdate: no -> deleting data from MSA tables", file = logfile)
    SQL <- paste("DELETE FROM sequence_selected", 
                 "WHERE", wrapSQL(gene, "locus"))
    dbSendQuery(conn, SQL)
  } 
  
  ## Select accessions that have an Expect value of at least 11
  ## and a coverage of al least 20%
  ## Note: sstrand = minus: cannot be rc'd, they are just bad
  ## -----------------------------------------------------------
  taxa <- paste("SELECT taxon, acc, length, mismatch, evalue, qcovs, sstrand",
               "FROM sequence", 
               "WHERE", wrapSQL(gene, "locus", "="),
               "AND evalue <= 11",
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
  chooseAcc <- function(tab, spec){
    tab <- tab[tab$taxon == spec, ]
    # if ("plus" %in% tab$sstrand){
    #   tab <- tab[tab$sstrand == "plus", ] # do not use minus if plus is available
    # }
    tab <- tab[tab$evalue == min(tab$evalue), ] # lowest E-value and
    tab <- tab[tab$qcovs == max(tab$qcovs), ] # greatest coverage
    tab <- tab[tab$length == max(tab$length), ] # longest alignment
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
  
  ## Write choosen accessions to MSA table
  ## -------------------------------------
  SQL <- paste("SELECT acc, taxon, sequence",
               "FROM sequence", 
               "WHERE", wrapSQL(selected, "acc", "=", "OR", by = 500))
  seqs <- lapply(SQL, dbGetQuery, conn = conn)
  seqs <- do.call(rbind, seqs)
  slog("\n", nrow(seqs), " sequences selected", sep = "", file = logfile)
  r <- range(sapply(seqs$sequence, nchar))
  slog("\nSequence lengths: ", r[1],"-", r[2], " bp", sep = "", file = logfile)
  
  ## Do reverse-complements 
  ## ----------------------
  # if ("minus" %in% chosen_acc$sstrand){
  #   id <- seqs$acc %in% chosen_acc$acc[chosen_acc$sstrand == "minus"]
  #   seqs$sequence[id] <- rcString(seqs$sequence[id])
  # }
  
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
  
  
  ## write to FASTA files
  ## ---------------------
  obj <- dbReadMSA(x)
  write.fas(obj, paste0("msa/", gene, ".fas"))
  
  dbDisconnect(conn)
  slog("\n\nSTEP F finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), "\n", file = logfile)
  dbProgress(x, "step_f", "success")
}
