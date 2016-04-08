## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-02-03)

compareToRef <- function(megProj, spec, reference){
  
  ## PARAMETERS
  ## ----------
  gene <- megProj@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  align.exe <- megProj@align.exe
  max.bp <- megProj@params@max.bp
  logfile <- paste(gene, "stepE.log", sep = "-")
  
  conn <- dbconnect(megProj)
  
  cat("\n", spec[1])
  
  obj <- dbReadDNA(megProj, acc.tab, spec[1], FALSE, max.bp)
  if ( is.list(obj) ) stop("sequences of ", spec[1], " not aligned")
  
  if ( length(spec) == 2 ){
    
    ## select appropriate reference sequences
    ## --------------------------------------
    br <- spec[2]
    ## for small groups there might not be a reference
    ## then choose the least distant reference available
    if ( !br %in% rownames(reference) ){
      obj <- c(as.list(reference), as.list(obj))
      obj <- mafft(obj, method = "auto", path = align.exe)
      d <- dist.dna(obj, model = "raw", pairwise.deletion = TRUE, 
                    as.matrix = TRUE)
      thisMean <- function(d, bn, sn){
        mean(d[bn, grep(sn, colnames(d))], na.rm = TRUE)
      }
      br <- sapply(names(reference), thisMean, d = d, sn = spec[1])
      br <- names(br)[which.min(br)]
      reference <- reference[br]
    }
  } else {
    br  <- "^REF_"
  }
  
  # profile alignment is far quicker
  obj <- mafft.merge(list(reference, obj), 
                     mafft.exe = megProj@align.exe,
                     thread = 1)
  
  ## handle reverse complements
  ## --------------------------
  rc <- grep("^_R_", rownames(obj))
  if ( length(rc) > 0 ){
    spec.id <- grep(spec[1], rownames(obj))
    rownames(obj) <- gsub("^_R_", "", rownames(obj))
    ## realignment if only subset of sequences concerned
    if ( length(spec.id) > length(rc) ){
      gt <- splitGiTaxon(rownames(obj)[rc])
      rc <- paste("UPDATE", acc.tab,
                  "SET status = 'aligned-RC'",
                  "WHERE", sql.wrap(gt[, 1], term = "gi", BOOL = NULL),
                  "AND", sql.wrap(gt[, 2], term = "taxon", BOOL = NULL))
      dbSendQuery(conn, rc)
      rc <- mafft(obj[spec.id, ], path = megProj@align.exe)
      dbWriteDNA(conn, acc.tab, rc)
      
      ## all sequences concerned
    } else {
      rc <- obj[rc, ]
      rc <- deleteEmptyCells(rc, TRUE)
      dbWriteDNA(conn, acc.tab, rc, status = "aligned-RC")
    } 
  }
  
  ## crop sequences
  ## --------------
  n <- ncol(obj)
  obj <- cropToReference(obj)
  if ( ncol(obj) < n ){
    dbWriteDNA(conn, acc.tab, obj[-grep("REF", rownames(obj)), ],
               status = "aligned")
  }
  
  ## calculate distance matrix
  ## -------------------------
  d <- dist.dna(obj, model = "raw", 
                pairwise.deletion = TRUE, 
                as.matrix = TRUE)
  ## assign distance = 1 to seqs that do
  ## not oberlap with reference seq
  d[is.na(d)] <- 1 
  
  ## convert into identity matrix
  ## ----------------------------
  d <- 1 - d
  d <- d[, grep(br, colnames(d)), drop = FALSE]
  d <- d[-grep(br, rownames(d)), , drop = FALSE]
  closest.ref <- colnames(d)[apply(d, 1, which.max)]
  d <- apply(d, 1, max)
  cc <- data.frame(d, closest.ref, stringsAsFactors = FALSE)
  
  ## calculate coverage in relation to next similar
  ## reference sequence
  ## ------------------
  cr <- unique(cc$closest.ref)
  covNextSim <- function(r, cc){
    cr <- c(r, rownames(cc)[cc$closest.ref == r])
    cr <- deleteEmptyCells(obj[cr, ], quiet = TRUE)
    cr <- coverage(cr)
    cr[-1]/cr[1] # first element is reference
  }
  cv <- lapply(cr, covNextSim, cc)
  cv <- unlist(cv)
  
  ## write results to database
  ## -------------------------
  d <- data.frame(splitGiTaxon(names(d)), 
                  identity = d,
                  coverage = cv)
  SQL <- paste("UPDATE", acc.tab, 
               "SET", wrapSQL(d$identity, term = "identity", 
                              operator = "=", boolean = NULL), 
               ",", wrapSQL(d$coverage, term = "coverage", 
                            operator = "=", boolean = NULL),
               "WHERE", wrapSQL(d$gi, term = "gi", 
                                operator = "=", boolean = NULL), 
               "AND", wrapSQL(d$taxon, term = "taxon", 
                              operator = "=", boolean = NULL))
  lapply(SQL, dbSendQuery, conn = conn)
  dbDisconnect(conn)
} 