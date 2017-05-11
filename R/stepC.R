## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-02-22)

#' @export
#' @import DBI
#' @import snow
#' @import snowfall
#' @importFrom snow setDefaultClusterOptions

stepC <- function(x){
  
  start <- Sys.time()
  pipeline <- TRUE
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("undefined locus not allowed")
  
  ## check if previous step has been run
  ## -----------------------------------
  status <- dbProgress(x)
  if (status$step_b == "pending") {
    stop("the previous step has not been called yet")
  }
  if (status$step_b == "error") {
    stop("the previous step has terminated with an error")
  }
  if (status$step_b == "failure") {
    slog("\nNo data from upstream available - quitting", file = "")
    dbProgress(x, "step_c", "failure")
    return()
  }
  if (status$step_b == "success") {
    dbProgress(x, "step_c", "error")
  }
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tax <- dbReadTaxonomy(x)
  align.exe <- x@align.exe
  max.bp <- x@params@max.bp
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepC.log")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""), 
       "\n\nSTEP C: alignment of conspecific sequences\n",
       file = logfile)
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x)
  
  ## clear results from previous runs of steps D
  ## -------------------------------------------
  SQL <- paste("UPDATE", acc.tab, 
               "SET status = 'raw'",
               "WHERE status ~ 'reference'")
  dbSendQuery(conn, SQL) 
  
  ## get table of species available in database
  ## ------------------------------------------
  SQL <- paste("SELECT taxon AS spec,", 
               "count(taxon) AS n,",
               "min(char_length(dna)) = max(char_length(dna)) AS aligned",
               "FROM", acc.tab, 
               "WHERE status !~ 'excluded'",
               "GROUP BY taxon",
               "ORDER BY taxon") 
  tax <- dbGetQuery(conn, SQL)
  
  ## no sequences in table: inform and quit
  ## --------------------------------------
  if (!nrow(tax)) {
    dbProgress(x, "step_c", "failure")
    dbDisconnect(conn)
    slog("no sequences - try to rerun stepB", file = logfile)
    slog("\n\nSTEP C finished", file = logfile)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), file = logfile)
    return()
  }
  
  ## mark single-sequence species in <status>
  ## ----------------------------------------
  singles <- tax$spec[tax$n == 1]
  slog("\n", nrow(tax), "species in table", acc.tab,
       "\n", length(singles), "species have 1 accession",
       "\n", nrow(tax) - length(singles), "species have > 1 accession",
       file = logfile)
  SQL <- paste("UPDATE", acc.tab,
               "SET status = 'single'",
               "WHERE", wrapSQL(singles, "taxon", "=", NULL),
               "AND status !~ 'excluded'")
  lapply(SQL, dbSendQuery, conn = conn)
  
  ## mark species that have already been aligned
  ## -------------------------------------------
  aligned <- tax$spec[tax$n > 1 & tax$aligned]
  slog("\n", length(aligned), "species are already aligned",
       file = logfile)
  SQL <- paste("UPDATE", acc.tab,
               "SET status = 'aligned'",
               "WHERE", wrapSQL(aligned, "taxon", "=", NULL),
               "AND status !~ 'excluded'")
  lapply(SQL, dbSendQuery, conn = conn)
  
  ## select species that need to be aligned
  ## --------------------------------------
  spec <- tax$spec[tax$n > 1 & !tax$aligned]
  slog("\n", length(spec), "species need to be aligned\n\n", 
       file = logfile)
  
  ## aligning -- either sequential or parallel
  ## -----------------------------------------
  if ( length(spec) > 0 ) {
    cpus <- x@params@cpus
    if ( length(spec) < cpus | !x@params@parallel ){
      lapply(spec, alignSpecies, megProj = x)
    } else {
      sfInit(parallel = TRUE, cpus = cpus, 
                       type = x@params@cluster.type)
      sfLibrary("megaptera", character.only = TRUE)
      megProj <- x
      sfExport("spec", "megProj", "acc.tab", 
                         "max.bp", "align.exe", "logfile")
      sfLapply(x = spec, fun = alignSpecies, megProj = megProj)
      sfStop()
    }
  } 
  dbDisconnect(conn)
  
  slog("\n\nSTEP C finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), 
       "\n", file = logfile)
  dbProgress(x, "step_c", "success")
}
