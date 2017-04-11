## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-03-28)

#' @export

stepF <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("please define locus")
  
  ## check if previous step has been run
  ## -----------------------------------
  status <- dbProgress(x)
  if (status$step_e == "pending") {
    stop("the previous step has not been called yet")
  }
  if (status$step_e == "error") {
    stop("the previous step has terminated with an error")
  }
  if (status$step_e == "failure") {
    slog("\nNo data from upstream available - quitting", file = "")
    dbProgress(x, "step_f", "failure")
    return()
  }
  if (status$step_e == "success") {
    dbProgress(x, "step_f", "error")
  }
  
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <-x@taxon@tip.rank
  msa.tab <- paste(tip.rank, gsub("^_", "", gene), sep = "_")
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
  
  ## set threshold values
  ## --------------------
  min.identity <- x@locus@min.identity
  if (min.identity < 0){
    min.identity <- optimizeIdentity(x, quiet = TRUE)
    min.identity <- min.identity$suggested - 1e-06
  }
  slog("\n.. threshold of minimum identity:", min.identity, 
       file = logfile)
  min.coverage <- x@locus@min.coverage
  if (min.coverage < 0){
    min.coverage <- optimizeCoverage(x, quiet = TRUE)
    min.coverage <- min.coverage$suggested - 1e-06
  }
  slog("\n.. threshold of minimum coverage:", min.coverage, 
       file = logfile)
  
  ## open database connection
  conn <- dbconnect(x)
  
  ## TO DO: nachfolgenden Code sinnvoll einbetten!
  if (!x@update){
    dbRemoveTable(conn, msa.tab)
  } 
  
  ## UNDER DEVOLOPMENT: check if update necessary
  ## --------------------------------------------
  SQL <- paste("SELECT taxon, gi",
               "FROM", acc.tab, 
               "WHERE npos <=", max.bp,
               "AND identity >=", min.identity,
               "AND coverage >=", min.coverage,
               "AND status !~ 'excluded|too'",
               "ORDER BY (taxon, gi)")
  taxa <- dbGetQuery(conn, SQL)
  if (!nrow(taxa)) {
    dbDisconnect(conn)
    slog(paste("\n.. WARNING: no sequences comply with min.identity=", 
               min.identity, " and min.coverage=", min.coverage, sep = ""), 
         file = logfile)
    dbProgress(x, "step_f", "failure")
    return()
  }
  if (tip.rank == "genus"){
    taxa$taxon <- strip.spec(taxa$taxon)
  }
  taxa <- split(taxa$gi, f = taxa$taxon)
  n <- sapply(taxa, length)
  taxa <- lapply(taxa, md5)
  taxa <- do.call(rbind, taxa)
  colnames(taxa) <- "md5"
  taxa <- data.frame(taxon = rownames(taxa), n, taxa, stringsAsFactors = FALSE)
  slog(paste("\n.. ", nrow(taxa), " species found in table '", 
             acc.tab, "'", sep = ""), file = logfile)
  
  if (dbExistsTable(conn, msa.tab)){
    
    ## stepF was run before, now check 
    ## if table needs updating, i.e. ...
    in_tab <- paste("SELECT taxon, md5",
                    "FROM", msa.tab)
    in_tab <- dbGetQuery(conn, in_tab)
    
    ## are there old sequences?
    old.spec <- !in_tab$taxon %in% taxa$taxon
    if (any(old.spec)){
      delete <- in_tab$taxon[old.spec]
      update_taxa <- paste("\n  -", delete)
      update_taxa <- paste(update_taxa, collapse = "")
      slog("\n.. removing", length(delete), "species:", 
           update_taxa, file = logfile)
      
      delete <- c(
        paste("DELETE FROM", msa.tab, 
              "WHERE", wrapSQL(delete, 
                               operator = "=",
                               boolean = NULL)),
        paste("UPDATE", msa.tab,
              "SET status = 'raw'",
              "WHERE status !~ 'excluded'"))
      lapply(delete, dbSendQuery, conn = conn)
    }
    ## ... are there new species?
    new.spec <- !taxa$spec %in% in_tab$taxon
    ## .. are there new sequences?
    new.seq <- !taxa$md5 %in% in_tab$md5
    
    if (any(new.spec) | any(new.seq)){
      
      ## table needs updating
      ## new species
      if (any(new.spec)){
        spec.id <- which(new.spec)
        update_taxa <- paste("\n  -", taxa[spec.id, "taxon"])
        update_taxa <- paste(update_taxa, collapse = "")
        slog("\n.. adding", length(spec.id), "new species:", 
             update_taxa, file = logfile)
      } else {
        spec.id <- NULL
      }
      ## new sequences
      if (any(new.seq)){
        seq.id <- which(new.seq)
        update_taxa <- paste("\n  -", taxa[seq.id, "taxon"])
        update_taxa <- paste(update_taxa, collapse = "")
        slog(paste0("\n.. set of available sequences has changed for ", 
                   length(seq.id), " ", tip.rank, ": ", update_taxa), 
             file = logfile)
      } else {
        seq.id <- NULL
      }
      ## subset of species that needs updating
      taxa <- taxa[union(spec.id, seq.id), ]
      
      ## reset status for species to update
      SQL <- paste("UPDATE", acc.tab,
                   "SET status = 'aligned'",
                   "WHERE (status = 'too short (from reference)'",
                   "OR status = 'too distant (from reference)')",
                   "AND", wrapSQL(taxa$spec, term = "taxon",
                                  operator = "=", 
                                  boolean = NULL))
      lapply(SQL, dbSendQuery, conn = conn)
      
      ## remove alignments files from previous runs
      fn <- list.files(pattern = paste0(gene, ".+[fas|phy|nex]$"))
      file.remove(fn)
      
      ## update <status> for accessions *NOT* selected
      ## ---------------------------------------------
      SQL <- c(paste("UPDATE", acc.tab, 
                     "SET status='excluded (too long)'", ## see stepC!
                     "WHERE npos >", max.bp),
               paste("UPDATE", acc.tab, 
                     "SET status='too distant (from reference)'",
                     "WHERE identity <", min.identity),
               paste("UPDATE", acc.tab, 
                     "SET status='too short (from reference)'",
                     "WHERE coverage <", min.coverage))
      lapply(SQL, dbSendQuery, conn = conn)
      
    } else {
      
      ## no updating necessary
      ## ---------------------
      dbDisconnect(conn)
      slog(paste("\n.. table '", msa.tab, 
                 "' is up to date - nothing to do\n", sep = ""), 
           file = logfile)
      dbProgress(x, "step_f", "success")
      return()
    }
  } else {
    
    slog("\n.. running step F for the first time", file = logfile)
    ## create table
    slog(paste("\n.. creating table '", msa.tab, "' ..", sep = ""), 
         file = logfile)
    SQL <- paste0(msa.tab, "_pk")
    SQL <- paste("CREATE TABLE", msa.tab, 
                 "(taxon text NOT NULL,",
                 "n integer,",
                 "md5 character(32),",
                 "subtree text,",
                 "status text,",
                 "npos integer NOT NULL,", 
                 "dna text NOT NULL,",
                 "masked text,",
                 "CONSTRAINT", SQL, "PRIMARY KEY (taxon))")
    dbSendQuery(conn, SQL)
  }
  
  slog("\n.. calculate majority-rule consensus sequences ..", 
       file = logfile)
  
  ## consensus -- either sequential or parallel
  ## ------------------------------------------
  if (nrow(taxa)) {
    cpus <- x@params@cpus
    if (nrow(taxa) < cpus | !x@params@parallel){
      apply(taxa, 1, speciesConsensus, megProj = x)
    } else {
      sfInit(parallel = TRUE, cpus = cpus, 
             type = x@params@cluster.type)
      sfLibrary("megaptera", character.only = TRUE) 
      megProj <- x
      sfExport("megProj", "taxa")
      sfApply(taxa, speciesConsensus, margin = 1, megProj = megProj)
      sfStop()
    }
  }
  
  ## write to FASTA files
  ## ---------------------
  obj <- dbReadDNA(x, msa.tab)
  write.fas(obj, paste0("msa/", gene, ".fas"))
  
  dbDisconnect(conn)
  slog("\n\nSTEP F finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), "\n", file = logfile)
  dbProgress(x, "step_f", "success")
}
