## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-04-08)

stepD <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  if ( x@locus@kind == "undefined" ) stop("undefined locus not allowed")
  
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  align.exe <- x@align.exe
  max.bp <- x@params@max.bp
  reference.max.dist <- x@params@reference.max.dist
  min.seqs.reference <- x@params@min.seqs.reference
  
  ## iniate logfile
  ## --------------
  logfile <- paste(gene, "stepD.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP D: reference sequence\n", file = logfile)
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x)
  
  ## check if stepC has been run
  ## ---------------------------
  status <- paste("SELECT DISTINCT status",
                  "FROM", acc.tab)
  status <- dbGetQuery(conn, status)
  if ( "raw" %in% status$status ){
    dbDisconnect(conn)
    stop("stepC has not been run")
  }
  
  
  
  
  ## delete previous entries
  ## -----------------------
  if ( dbExistsTable(conn, "reference") ){
    SQL <- paste("DELETE FROM reference",
                 "WHERE", sql.wrap(gene, term = 'gene'))
    dbSendQuery(conn, SQL)
  }
  
  ## CASE 1: user-defined reference sequences
  ## ----------------------------------------
  if ( inherits(x@locus, "locusRef") ){
  
    slog("\n.. user-defined reference sequences ..", 
         file = logfile)
    ref <- x@locus@reference
    ref <- as.list(ref)
    if ( length(ref) > 1 ){
      ref <- mafft(ref, path = x@align.exe) # realign in any case!
      ref <- removeIsolatedTails(ref)
      ref <- trimEnds(ref, .75 * nrow(ref)) ## try to make coverage comparable
    } else {
      ref <- as.matrix(ref)
    }
    
    dbUpdateReference(conn, gene, ref)
    
  } else {
    
    ## CASE 2: pipeline-defined reference sequences
    ## --------------------------------------------
    slog("\n.. using heuristic algorithm", 
         "to create reference sequences ..", file = logfile)
    slog("\n.. minimum number of species required for reference calculation:", 
         min.seqs.reference, file = logfile)
    
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
    slog("\n.. calculate reference(s) for (sub-)groups at rank:", 
         rr, file = logfile)
    
    ## reset database if updating
    ## --------------------------
    if ( x@update & dbExistsTable(conn, "reference") ) {
      SQL <- c(paste("UPDATE", acc.tab, 
                     "SET status='aligned' WHERE status~'reference'"),
               paste("DELETE FROM reference WHERE", 
                     sql.wrap(gene, term = "gene")))
      lapply(SQL, dbSendQuery, conn = conn)
    }
    
    ## species, for which we have more than 1 acc.
    ## -------------------------------------------
    tax <- paste("SELECT taxon AS spec, count(taxon),", rr,
                 "FROM", acc.tab, 
                 "JOIN taxonomy ON (taxon = spec)", 
                 "WHERE status ~ 'aligned'",
                 "GROUP BY taxon,", rr,
                 "ORDER by taxon")
    tax <- dbGetQuery(conn, tax)
    
    if ( nrow(tax) == 0 ) {
      slog("\nWARNING: no species with > 1 accession available", file = logfile)
      dbDisconnect(conn)
      slog("\n\nSTEP D finished", file = logfile)
      td <- Sys.time() - start
      slog(" after", round(td, 2), attr(td, "units"), file = logfile)
      return(x)
    } else {
      slog("\n..", nrow(tax), "species have > 1 sequence", file = logfile)
    }
    
    ## reference.clades: updating or not ...
    ## ------------------------------------
    reference.clades <- unique(tax[[rr]])
    if ( dbExistsTable(conn, "reference") ){
      already <- names(dbReadReference(x))
      if ( !x@update & all(reference.clades %in% already) ){
        slog("\n.. reference sequences already calculated", 
             file = logfile)
        dbDisconnect(conn)
        return(x)
      } else {
        reference.clades <- setdiff(reference.clades, already)
      }
    } else {
      already <- NULL
    }
    
    slog("\n.. calculating reference sequences for", 
         length(reference.clades), "clades", 
         file = logfile)
    
    ref <- vector(mode = "list", length = length(reference.clades))
    names(ref) <- reference.clades
    for ( i in reference.clades ) {
      
      slog("\n -", i, file = logfile)
      
      spec <- tax[tax[[rr]] == i, "spec"]
      ali <- lapply(spec, dbReadDNA, tab.name = acc.tab, x = x)
      names(ali) <- gsub(" ", "_", spec)
      slog("\n  > number of species with > 1 accession:", 
           length(ali), file = logfile)
      
      ## identify species, which have (nearly) identical accessions
      ## ----------------------------------------------------------
      id <- data.frame(sapply(ali, getMaxDist))
      id <- split(rownames(id), f = id[, 1])
      d <- as.numeric(names(id))
      for ( j in d ){
        bt <- sort(unlist(id[d <= j]))
        if ( length(bt) >= min.seqs.reference ) break
      }
      ali <- ali[bt]
      bp <- ifelse(j == 0, "bp", "bp or less")
      slog("\n  > found", length(bt), "species with a distance of", j, bp,  
           file = logfile)
      
      ## check if there are species that have not been used 
      ## to create reference. Only create reference if new species are available
      ## -----------------------------------------------------------------------
      sql <- sql.wrap(names(ali), term = "taxon")
      sql <- paste("(", sql, ")", sep = "")
      sql <- paste("status ~ 'reference' AND", sql)
      SQL <- paste("SELECT taxon, status FROM", acc.tab, "WHERE", sql)
      SQL <- unique(dbGetQuery(conn, SQL))
      
      if ( length(grep("reference", SQL$status)) == length(ali) & 
             i %in% already ){ # it is possible that status column contains
        # 'reference', but reference table has been
        # deleted since
        slog("\n- reference has already been calculated", file = logfile)
        next
      } 
      
      ## sequences used for reference calculation
      u <- unlist(lapply(ali, rownames))
      u <- splitGiTaxon(u)
      
      ## ... and make a consensus sequence for each of them.
      ## ---------------------------------------------------
      ali <- lapply(ali, specCons)
      class(ali) <- "DNAbin"
      
      ## Then align 
      ## ----------
      if ( length(ali) > 1 ) {
        ali <- mafft(ali, method = "auto", path = align.exe)
      } else {
        ali <- as.matrix(ali)  
      }
      cv <- coverage(ali)
      ali <- ali[names(cv)[cv >= .5], ]
      ali <- trimEnds(ali, min.n.seq = 1)
      
      ## check distance matrix
      ## ---------------------
      if ( nrow(ali) > 1 ){
        slog("\n  > check genetic distances ...", file = logfile)
        d <- dist.dna(ali, model = "raw", as.matrix = TRUE,
                      pairwise.deletion = TRUE)
        diag(d) <- NA
        e <- apply(d, 1, min, na.rm = TRUE)
        id <- names(e)[e < reference.max.dist]
        if ( length(id) > 0 ){
          slog("\n   ", length(id), "sequences with maximum distance of",
               reference.max.dist, file = logfile)
          ali <- ali[names(e)[e < reference.max.dist], ]
          ali <- deleteEmptyCells(ali, quiet = TRUE)
        } else {
          slog("\n", length(id), "sequences conform to maximum reference distance of",
               reference.max.dist, "\n CAUTION: reference distance will be ignored", file = logfile)
          id <- rownames(ali)
        }
      } else {
        id <- rownames(ali) ## only one species
      }
      
      
      ## add 'reference' to status
      ## -------------------------
      id <- u[u[, "taxon"] %in% rownames(ali), "gi"] 
      SQL <- paste("UPDATE", acc.tab, 
                   "SET status='aligned-reference'",
                   "WHERE", sql.wrap(id, term = "gi"))
      dbSendQuery(conn, SQL)
      
      ## create 'reference' sequence
      ## ---------------------------
      slog("\n  > creating reference sequence ...", file = logfile)
      ref[[i]] <- specCons(ali)
    } # end of FOR-loop over i
    
    ## Align reference sequences ... 
    ## -----------------------------
    class(ref) <- "DNAbin"
    if ( length(ref) > 1 ){
      ref <- mafft(ref, path = align.exe)
    } else {
      ref <- as.matrix(ref)
    }
    
    ## ... and insert into table or update table
    ## -----------------------------------------
    dbUpdateReference(conn, gene, ref)
  }
  dbDisconnect(conn)
  
  slog("\n\nSTEP D finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), 
       "\n", file = logfile)
  invisible(x)
}