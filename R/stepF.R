## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-09-19)

stepF <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  if ( x@locus@kind == "undefined" ) stop("please define locus")
  STATUS <- checkStatus(x)
  if ( !all(STATUS[1:5]) ){
    stop("step", names(STATUS)[min(which(!STATUS))] ,
         " has not been run")
  }
  
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- x@taxon@tip.rank
  msa.tab <- paste(tip.rank, gsub("^_", "", gene), sep = "_")
  align.exe <- x@align.exe
  max.bp <- x@params@max.bp
  
  ## iniate logfile
  ## --------------
  logfile <- paste(gene, "stepF.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  r <- ifelse(tip.rank == "gen", "genus", "species")
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP F: construct", r, "consensus sequences\n", 
       file = logfile)
  
  ## set threshold values
  ## --------------------
  min.identity <- x@locus@min.identity
  if ( min.identity < 0 ){
    min.identity <- optimizeIdentity(x, quiet = TRUE)
    min.identity <- min.identity$suggested - 1e-06
  }
  slog("\n.. threshold of minimum identity:", min.identity, 
       file = logfile)
  min.coverage <- x@locus@min.coverage
  if ( min.coverage < 0 ){
    min.coverage <- optimizeCoverage(x, quiet = TRUE)
    min.coverage <- min.coverage$suggested - 1e-06
  }
  slog("\n.. threshold of minimum coverage:", min.coverage, 
       file = logfile)
  
  ## open database connection
  conn <- dbconnect(x)
  
  ## UNDER DEVOLOPMENT: check if update necessary
  ## --------------------------------------------
  SQL <- paste("SELECT", tip.rank, "AS taxon, gi",
               "FROM taxonomy",
               "JOIN", acc.tab, "AS a", 
               "ON taxonomy.spec = a.taxon",
               "WHERE npos <=", max.bp,
               "AND identity >=", min.identity,
               "AND coverage >=", min.coverage,
               "AND status !~ 'excluded|too'",
               "ORDER BY (taxon, gi)")
  taxa <- dbGetQuery(conn, SQL)
  if ( nrow(taxa) == 0 ) {
    dbDisconnect(conn)
    slog(paste("\n.. WARNING: no sequences comply with min.identity=", 
               min.identity, " and min.coverage=", min.coverage, sep = ""), 
         file = logfile)
    return()
  }
  taxa <- split(taxa$gi, f = taxa$taxon)
  n <- sapply(taxa, length)
  taxa <- lapply(taxa, md5)
  taxa <- do.call(rbind, taxa)
  colnames(taxa) <- "md5"
  taxa <- data.frame(taxon = rownames(taxa), n, taxa)
  slog(paste("\n.. ", nrow(taxa), " species found in table '", 
             acc.tab, "'", sep = ""), file = logfile)
  
  if ( dbExistsTable(conn, msa.tab) ){
    
    ## stepF was run before, now check 
    ## if table needs updating, i.e. ...
    in.tab <- paste("SELECT", tip.rank, ", md5",
                    "FROM", msa.tab)
    in.tab <- dbGetQuery(conn, in.tab)
    
    ## are there old sequences?
    old.spec <- !in.tab$spec %in% taxa$taxon
    if ( any(old.spec) ){
      delete <- in.tab$spec[old.spec]
      update.spec <- paste("\n  -", delete)
      update.spec <- paste(update.spec, collapse = "")
      slog("\n.. removing", length(delete), "species:", 
           update.spec, file = logfile)
      
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
    new.spec <- !taxa$spec %in% in.tab$spec
    ## .. are there new sequences?
    new.seq <- !taxa$md5 %in% in.tab$md5
    
    if ( any(new.spec) | any(new.seq) ){
      
      ## table needs updating
      ## new species
      if ( any(new.spec) ){
        spec.id <- which(new.spec)
        update.spec <- paste("\n  -", taxa$spec[spec.id])
        update.spec <- paste(update.spec, collapse = "")
        slog("\n.. adding", length(spec.id), "new species:", 
             update.spec, file = logfile)
      } else {
        spec.id <- NULL
      }
      ## new sequences
      if ( any(new.seq) ){
        seq.id <- which(new.seq)
        update.spec <- paste("\n  -", taxa$spec[seq.id])
        update.spec <- paste(update.spec, collapse = "")
        slog(paste("\n.. set of available sequences has changed for ", 
                   length(seq.id), " species:", update.spec,
                   sep = ""), file = logfile)
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
      fn <- list.files(pattern = paste(gene, ".+[phy|nex]$", sep = ""))
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
      return()
    }
  } else {
    
    slog("\n.. running step F for the first time", file = logfile)
    ## create table
    slog(paste("\n.. creating table '", msa.tab, "' ..", sep = ""), 
         file = logfile)
    SQL <- paste(msa.tab, "_pk", sep = "")
    SQL <- paste("CREATE TABLE", msa.tab, 
                 "(", tip.rank, "text NOT NULL,",
                 "n integer,",
                 "md5 text,",
                 "subtree text,",
                 "status text,",
                 "npos integer NOT NULL,", 
                 "dna text NOT NULL,",
                 "masked text,",
                 "CONSTRAINT", SQL, "PRIMARY KEY (", tip.rank, "))")
    dbSendQuery(conn, SQL)
  }
  
  slog("\n.. calculate majority-rule consensus sequences ..", 
       file = logfile)
  
  ## consensus -- either sequential or parallel
  ## ------------------------------------------
  if ( nrow(taxa) > 0 ) {
    cpus <- x@params@cpus
    if ( nrow(taxa) < cpus | !x@params@parallel ){
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
  
  dbDisconnect(conn)
  slog("\n\nSTEP F finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
  invisible(x)
}
