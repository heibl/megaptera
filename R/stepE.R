## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-04-08)

stepE <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  if ( x@locus@kind == "undefined" ) stop("undefined locus not allowed")
  
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  
  ## iniate logfile
  ## --------------
  logfile <- paste(gene, "stepE.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP E: calculate genetic distances from reference", 
       file = logfile)
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
  
  ## check if stepC has been run
  ## ---------------------------
  status <- paste("SELECT DISTINCT status",
                  "FROM", acc.tab)
  status <- dbGetQuery(conn, status)
  if ( "raw" %in% status$status ){
    dbDisconnect(conn)
    stop("stepC has not been run")
  }
  
  ## read reference sequences
  ## ------------------------
  reference <- dbReadReference(x)
  if ( is.logical(reference) ) {
    slog("\nWARNING: no reference sequences available", file = logfile)
    dbDisconnect(conn)
    slog("\n\nSTEP E finished", file = logfile)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), file = logfile)
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
  if ( !is.matrix(reference) ){
    reference <- mafft(reference, path = x@align.exe)
  } 
  
  ## CASE 1: user-defined reference sequences
  ## ----------------------------------------
  if ( inherits(x@locus, "locusRef") ){
    
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
    
    ## set taxonomic rank for 
    ## reference sequence calculation
    ## ------------------------------
    rr <- x@taxon@reference.rank
    if ( rr == "auto" ){
      tax <- dbReadTaxonomy(x)
      species.list <- unique(is.Linnean(unlist(x@taxon@ingroup)))
      if ( length(species.list) > 1 ) stop("names of species and higher taxa must not be mixed")
      if ( species.list ){
        rr <- apply(tax, 2, function(x) length(unique(x)))
        rr <- head(rr, -1) ## delete synonyms column
        rr <- names(rr)[max(which(rr == 1))]
      } else {
        rr <- apply(tax, 2, grep, pattern = paste("^", "$", 
                                                  sep = x@taxon@ingroup))
        rr <- names(rr)[sapply(rr, length) > 0]
      }
    }
    
    ## table of species names + reference sequence
    ## -------------------------------------------
    slog("\n.. reading species names: ", file = logfile)
    tax <- paste("SELECT DISTINCT spec,", rr, "AS ref",
                 "FROM", acc.tab, "JOIN taxonomy ON spec = taxon ",
                 "WHERE status !~ 'excluded|too'", 
                 "AND identity IS NULL",
                 "AND npos <=", x@params@max.bp,
                 "ORDER BY spec")
    tax <- dbGetQuery(conn, tax)
    slog(nrow(tax), "found ...", file = logfile)
    
    if ( nrow(tax) == 0 ) {
      slog("\n.. database is up to date -- nothing to do", file = logfile)
      dbDisconnect(conn)
      #     x$evaluate <- FALSE
      slog("\n\nSTEP E finished", file = logfile)
      td <- Sys.time() - start
      slog(" after", round(td, 2), attr(td, "units"), file = logfile)
      return(x)
    }
    
    ## Some reference taxa might not have their own reference,
    ## because they consist of only 1 single sequence.
    ## Here we use a hack ignoring that we do not know, which 
    ## reference is closest if there are more than 1 reference
    ## -------------------------------------------------------
    miss.ref <- which(!tax$ref %in% rownames(reference))
    if ( length(miss.ref) > 0 ){
      tax$ref[miss.ref] <- rownames(reference)[1]
    }
    
    ## convert data frame to list
    ## -------------------------
    tax <- paste(tax$spec, tax$ref, sep = "Zzz")
    tax <- strsplit(tax, "Zzz")
    
    ## select the 'best' sequences
    ## ---------------------------
    slog("\n.. calculating genetic distance from 'reference' for ...", 
         file = logfile)
    
    ## distance from reference -- either sequential or parallel
    ## --------------------------------------------------------
    if ( length(tax) > 0 ) {
      cpus <- x@params@cpus
      if ( length(tax) < cpus | !x@params@parallel ){
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
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
  invisible(x)
}
