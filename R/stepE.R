## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-04-11)

#' @export

stepE <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("undefined locus not allowed")
  
  ## check if previous step has been run
  ## -----------------------------------
  status <- dbProgress(x)
  if (status$step_d == "pending") {
    stop("the previous step has not been called yet")
  }
  if (status$step_d == "error") {
    stop("the previous step has terminated with an error")
  }
  if (status$step_d == "failure") {
    slog("\nNo data from upstream available - quitting", file = "")
    dbProgress(x, "step_e", "failure")
    return()
  }
  if (status$step_d == "success") {
    dbProgress(x, "step_e", "error")
  }
  
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- match.arg(x@taxon@tip.rank, c("species", "genus"))
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepE.log")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP E: calculate genetic distances from reference", 
       file = logfile)

  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
  
  ## read reference sequences
  ## ------------------------
  reference <- dbReadReference(x)
  if (is.logical(reference)) {
    slog("\nWARNING: no reference sequences available", file = logfile)
    dbDisconnect(conn)
    slog("\n\nSTEP E finished", file = logfile)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), file = logfile)
    dbProgress(x, "step_e", "failure")
    return()
  }
  
  ## clear reference column if all entries are to be updated
  ## -------------------------------------------------------
  ##  TO DO: replace by forced update/tabula rasa
#   if ( x@update ) {
#     dbSendQuery(conn, paste("UPDATE", acc.tab, 
#                             "SET identity = NULL,", 
#                             "coverage = NULL"))
#     SQL <- paste("UPDATE", acc.tab,
#                  "SET status = 'aligned'",
#                  "WHERE status = 'too short (from reference)'",
#                  "OR status = 'too distant (from reference)'")
#     dbSendQuery(conn, SQL)
#   }
  
  ## align references (if necessary)
  ## -------------------------------
  if (!is.matrix(reference)){
    reference <- mafft(reference, exec = x@align.exe)
  } 
  
  ## CASE 1: user-defined reference sequences
  ## ----------------------------------------
  if (inherits(x@locus, "locusRef")){
    
    slog("\n.. user-defined reference sequences ..", 
         file = logfile)
    
    ## mark reference sequences
    ## ------------------------
    rownames(reference) <- strip.infraspec(rownames(reference))
    rownames(reference) <- paste("REF", rownames(reference), sep = "_")
    
    ## table of species names + reference sequence
    ## -------------------------------------------
    slog("\n.. reading species names: ", file = logfile)
    tax <- paste("SELECT DISTINCT taxon",
                 "FROM", acc.tab, 
                 "WHERE status !~ 'excluded|too'",
                 "AND npos <=", x@params@max.bp,
                 "AND identity IS NULL",
                 "ORDER BY taxon")
    tax <- dbGetQuery(conn, tax)
    slog(nrow(tax), "found ...", file = logfile)
    
    ## distance from reference -- either sequential or parallel
    ## --------------------------------------------------------
    if ( nrow(tax) > 0 ) {
      cpus <- x@params@cpus
      if ( nrow(tax) < cpus | !x@params@parallel ){
        lapply(tax$taxon, compareToRef, 
              megProj = x, reference = reference)
      } else {
        slog("\n", file = logfile)
        sfInit(parallel = TRUE, cpus = cpus, 
               type = x@params@cluster.type)
        sfLibrary("megaptera", character.only = TRUE)
        sfExport("reference", "x")
        sfLapply(tax$taxon, compareToRef, megProj = x, reference = reference)
        sfStop()
      }
    }
  } else {
    
    ## CASE 2: pipeline-defined reference sequences
    ## --------------------------------------------
    slog("\n.. pipeline-defined reference sequences ..", 
         file = logfile)
    
    ## get a list of species to work on ...
    ## ------------------------------------
    slog("\n.. reading species names: ", file = logfile)
    tax <- dbReadTaxonomy(x)
    tax <- paste("SELECT DISTINCT taxon",
                 "FROM", acc.tab, 
                 "WHERE status !~ 'excluded|too'", 
                 "AND identity IS NULL",
                 "AND npos <=", x@params@max.bp,
                 "ORDER BY taxon")
    tax <- dbGetQuery(conn, tax)
    if (!nrow(tax)) {
      slog("\n.. database is up to date -- nothing to do", file = logfile)
      dbDisconnect(conn)
      #     x$evaluate <- FALSE
      slog("\n\nSTEP E finished", file = logfile)
      td <- Sys.time() - start
      slog(" after", round(td, 2), attr(td, "units"), file = logfile)
      dbProgress(x, "step_e", "success")
      return(x)
    }
    
    ## ... and find their reference clade
    ## ---------------------------------
    rr <- x@taxon@reference.rank
    if (rr == "auto"){
      gt <- comprehensiveGuidetree(x, tip.rank = "species", subset = tax$taxon)
      refc <- gt$edge[gt$edge[, 1] == (Ntip(gt) + 1), 2]
      refc <- lapply(refc, descendants, phy = gt, labels = TRUE)
      names(refc) <- sapply(refc, taxdumpMRCA, x = x, tip.rank = "species")
      tax <- data.frame(ref = rep(names(refc), sapply(refc, length)), 
                        taxon = gsub("_", " ", unlist(refc)),
                        stringsAsFactors = FALSE)
      refc <- names(refc)
    } else {
      refc <- dbReadTaxonomy(x)
      tax <- data.frame(ref = taxdumpHigherRank(refc, tax$taxon, rr),
                        tax)
      refc <- unique(tax$rr)
    }
    slog(nrow(tax), "found ...", file = logfile)
  
    ## Some reference taxa might not have their own reference,
    ## because they consist of only 1 single sequence.
    ## Here we use a hack ignoring that we do not know, which 
    ## reference is closest if there are more than 1 reference
    ## -------------------------------------------------------
    miss.ref <- which(!tax$ref %in% rownames(reference))
    if (length(miss.ref)){
      tax$ref[miss.ref] <- rownames(reference)[1]
    }
    
    ## convert data frame to list
    ## -------------------------
    tax <- paste(tax$taxon, tax$ref, sep = "Zzz")
    tax <- strsplit(tax, "Zzz")
    
    ## select the 'best' sequences
    ## ---------------------------
    slog("\n.. calculating genetic distance from 'reference' for ...", 
         file = logfile)
    
    ## distance from reference -- either sequential or parallel
    ## --------------------------------------------------------
    if (length(tax)) {
      cpus <- x@params@cpus
      if (length(tax) < cpus | !x@params@parallel){
        lapply(tax, compareToRef, 
              megProj = x, reference = reference)
      } else {
        slog("\n", file = logfile)
        sfInit(parallel = TRUE, cpus = cpus, 
               type = x@params@cluster.type)
        sfLibrary("megaptera", character.only = TRUE)
        sfExport("reference", "x")
        sfLapply(tax, compareToRef, megProj = x, reference = reference)
        sfStop()
      }
    }
  }
  dbDisconnect(conn)
  slog("\n\nSTEP E finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), "\n", file = logfile)
  dbProgress(x, "step_e", "success")
}
