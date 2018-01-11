## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-11-28)

#' @export

stepD <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("undefined locus not allowed")
  
  ## check if previous step has been run
  ## -----------------------------------
  status <- dbProgress(x)
  if (status$step_c == "pending") {
    stop("the previous step has not been called yet")
  }
  if (status$step_c == "error") {
    stop("the previous step has not been called yet")
  }
  if (status$step_c == "failure") {
    slog("\nNo data from upstream available - quitting", file = "")
    dbProgress(x, "step_d", "failure")
    return()
  }
  if (status$step_c == "success") {
    dbProgress(x, "step_d", "error")
  }
  
  
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- match.arg(x@taxon@tip.rank, c("species", "genus"))
  align.exe <- x@align.exe
  max.bp <- x@params@max.bp
  reference.max.dist <- x@params@reference.max.dist
  min.seqs.reference <- x@params@min.seqs.reference
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepD.log")
  if (file.exists(logfile)) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP D: reference sequence\n", file = logfile, megProj = x)
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x)
  
  ## delete previous entries
  ## -----------------------
  if (dbExistsTable(conn, "reference")){
    SQL <- paste("DELETE FROM reference",
                 "WHERE", wrapSQL(gene, "gene", "="))
    dbSendQuery(conn, SQL)
  }
  
  ## CASE 1: user-defined reference sequences
  ## ----------------------------------------
  if (inherits(x@locus, "locusRef")){
    
    slog("\n.. user-defined reference sequences ..", 
         file = logfile)
    ref <- x@locus@reference
    ref <- as.list(ref)
    if ( length(ref) > 1 ){
      ref <- mafft(ref, exec = x@align.exe) # realign in any case!
      ref <- trimEnds(ref, .75 * nrow(ref)) ## try to make coverage comparable
    } else {
      ref <- as.matrix(ref)
    }
    
    dbUpdateReference(conn, gene, ref)
    
  } else {
    
    ## CASE 2: pipeline-defined reference sequences
    ## --------------------------------------------
    slog("\nUsing heuristic algorithm to create reference sequences", 
         file = logfile, megProj = x)
    slog("\nMinimum number of species required for reference calculation:", 
         min.seqs.reference, file = logfile, megProj = x)
    
    ## reset database if updating
    ## --------------------------
    if (x@update & dbExistsTable(conn, "reference")) {
      SQL <- c(paste("UPDATE", acc.tab, 
                     "SET status='aligned' WHERE status~'reference'"),
               paste("DELETE FROM reference WHERE", 
                     wrapSQL(gene, "gene", "=")))
      lapply(SQL, dbSendQuery, conn = conn)
    }
    
    ## do we have any ingroup taxa?
    ## ----------------------------
    tax <- paste("SELECT taxon, count(taxon) AS n",
                 "FROM", acc.tab,
                 "WHERE status !~ 'excluded|too'",
                 "AND", wrapSQL(x@params@max.bp, "npos", "<="),
                 "GROUP BY taxon")
    tax <- dbGetQuery(conn, tax)
    n_ingroup <- length(which(is.ingroup(x, tax$taxon)))
    if (n_ingroup < 3){
      slog("\nWARNING: less than 3 ingroup species available", file = logfile)
      dbDisconnect(conn)
      slog("\n\nSTEP D finished", file = logfile)
      td <- Sys.time() - start
      slog(" after", round(td, 2), attr(td, "units"), file = logfile)
      dbProgress(x, "step_d", "failure")
      return()
    } 
    
    ## which species have more than 1 accession?
    ## -----------------------------------------
    id <- tax$n > 1
    if (!nrow(tax[id, ])) {
      slog("\nWARNING: no species with > 1 accession available", file = logfile)
      dbDisconnect(conn)
      slog("\n\nSTEP D finished", file = logfile)
      td <- Sys.time() - start
      slog(" after", round(td, 2), attr(td, "units"), file = logfile)
      dbProgress(x, "step_d", "failure")
      return()
    } else {
      slog("\n..", length(id), "species have > 1 sequence", file = logfile)
    }
    
    ## determine reference clades (refc)
    ## ---------------------------------
    rr <- x@taxon@reference.rank
    if (rr == "auto"){
      gt <- comprehensiveGuidetree(x, tip.rank = "species", subset = tax$taxon)
      ref <- gt$edge[gt$edge[, 1] == (Ntip(gt) + 1), 2]
      ref <- lapply(ref, descendants, phy = gt, labels = TRUE)
      ref <- lapply(ref, gsub, pattern = "_", replace = " ")
      names(ref) <- sapply(ref, taxdumpMRCA, x = x, tip.rank = "species")
      ref <- data.frame(ref = rep(names(ref), sapply(ref, length)), 
                         taxon = gsub("_", " ", unlist(ref)), stringsAsFactors = FALSE)
      tax <- data.frame(ref, n = tax[match(ref$taxon, tax$taxon), "n"])
      refc <- unique(tax$ref)
    } else {
      refc <- dbReadTaxonomy(x)
      tax <- data.frame(ref = taxdumpHigherRank(refc, tax$taxon, rr),
                        tax)
      refc <- unique(tax$rr)
    }
    dbUpdateReference_Clade(x, tax) ## write to database
    
    ## select taxa with more than 1 accession
    ## --------------------------------------
    tax <- tax[tax$n > 1, ]
    refc <- intersect(tax$ref, refc)
  
    ## reference clades (refc): updating or not
    ## ----------------------------------------
    if (dbExistsTable(conn, "reference")){
      already <- names(dbReadReference(x))
      if (!x@update & all(refc %in% already)){
        slog("\n.. reference sequences already calculated", 
             file = logfile)
        dbDisconnect(conn)
        dbProgress(x, "step_d", "success")
        return()
      } else {
        refc <- setdiff(refc, already)
      }
    } else {
      already <- NULL
    }
    
    slog("\n.. calculating reference sequences for", 
         length(refc), "clades", file = logfile)
    ref <- vector(mode = "list", length = length(refc))
    names(ref) <- refc
    for (i in refc) {
      
      slog("\n -", i, file = logfile, megProj = x)
      
      spec <- tax[tax$ref == i, "taxon"]
      if (!length(spec)){
        slog("\nCAUTION: reference clade extended", 
             length(ali), file = logfile, megProj = x)
        spec <- tax[, "taxon"]
      }
      ali <- lapply(spec, dbReadDNA, x = x, tab.name = acc.tab, regex = FALSE)
      names(ali) <- gsub(" ", "_", spec)
      slog("\n  > number of species with > 1 accession:", 
           length(ali), file = logfile, megProj = x)
      
      ## identify species, which have (nearly) identical accessions
      ## ----------------------------------------------------------
      id <- data.frame(sapply(ali, getMaxDist))
      id <- split(rownames(id), f = id[, 1])
      d <- as.numeric(names(id))
      for (j in d){
        bt <- sort(unlist(id[d <= j]))
        if (length(bt) >= min.seqs.reference) break
      }
      ali <- ali[bt]
      bp <- ifelse(j == 0, "bp", "bp or less")
      slog("\n  > found", length(bt), "species with a distance of", j, bp,  
           file = logfile, megProj = x)
      
      ## check if there are species that have not been used 
      ## to create reference. Only create reference if new species are available
      ## -----------------------------------------------------------------------
      sql <- wrapSQL(names(ali), "taxon", "=")
      sql <- paste0("(", sql, ")")
      sql <- paste("status ~ 'reference' AND", sql)
      SQL <- paste("SELECT taxon, status FROM", acc.tab, "WHERE", sql)
      SQL <- unique(dbGetQuery(conn, SQL))
      
      if (length(grep("reference", SQL$status)) == length(ali) & 
           i %in% already){ # it is possible that status column contains
        # 'reference', but reference table has been
        # deleted since
        slog("\n- reference has already been calculated", file = logfile, megProj = x)
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
      if (length(ali) > 1) {
        ali <- mafft(ali, method = "auto", exec = align.exe)
      } else {
        ali <- as.matrix(ali)  
      }
      cv <- coverage(ali)
      ali <- ali[names(cv)[cv >= .5], ]
      ali <- trimEnds(ali, min.n.seq = 1)
      
      ## check distance matrix
      ## ---------------------
      if (nrow(ali) > 1 ){
        slog("\n  > check genetic distances ...", file = logfile, megProj = x)
        d <- dist.dna(ali, model = "raw", as.matrix = TRUE,
                      pairwise.deletion = TRUE)
        diag(d) <- NA
        e <- apply(d, 1, min, na.rm = TRUE)
        id <- names(e)[e < reference.max.dist]
        if ( length(id)){
          slog("\n   ", length(id), "sequences with maximum distance of",
               reference.max.dist, file = logfile, megProj = x)
          ali <- ali[names(e)[e < reference.max.dist], ]
          ali <- deleteEmptyCells(ali, quiet = TRUE)
        } else {
          slog("\n", length(id), "sequences conform to maximum reference distance of",
               reference.max.dist, "\n CAUTION: reference distance will be ignored", file = logfile, megProj = x)
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
                   "WHERE", wrapSQL(id, "gi", "="))
      dbSendQuery(conn, SQL)
      
      ## create 'reference' sequence
      ## ---------------------------
      slog("\n  > creating reference sequence ...", file = logfile, megProj = x)
      ref[[i]] <- specCons(ali)
    } # end of FOR-loop over i
    
    ## Align reference sequences ... 
    ## -----------------------------
    class(ref) <- "DNAbin"
    names(ref) <- gsub(" ", "_", names(ref))
    if (length(ref) > 1){
      ref <- mafft(ref, exec = align.exe)
    } else {
      ref <- as.matrix(ref)
    }
    
    ## ... and insert into table or update table
    ## -----------------------------------------
    dbUpdateReference(conn, gene, ref)
  }
  dbDisconnect(conn)
  
  slog("\n\nSTEP D finished", file = logfile, megProj = x)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), 
       "\n", file = logfile, megProj = x)
  dbProgress(x, "step_d", "success")
}