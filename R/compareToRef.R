## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-11-29)

#' @importFrom ips mafft.merge
#' @export

compareToRef <- function(megProj, spec, reference){
  
  ## PARAMETERS
  ## ----------
  gene <- megProj@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  align.exe <- megProj@align.exe
  max.bp <- megProj@params@max.bp
  logfile <- paste0("log/", gene, "-stepE.log")
  
  conn <- dbconnect(megProj)
  
  cat("\n", spec[1])
  
  obj <- dbReadDNA(megProj, tab.name = acc.tab, taxon = spec[1], 
                   regex = FALSE, max.bp = max.bp)
  if (is.list(obj)) stop("sequences of ", spec[1], " not aligned")
  
  if (length(spec) == 2){
    
    ## select appropriate reference sequences
    ## --------------------------------------
    br <- spec[2]
    ## for small groups there might not be a reference
    ## then choose the least distant reference available
    if ( !br %in% rownames(reference) ){
      obj <- c(as.list(reference), as.list(obj))
      obj <- mafft(obj, method = "auto", exec = align.exe)
      d <- dist.dna(obj, model = "raw", pairwise.deletion = TRUE, 
                    as.matrix = TRUE)
      thisMean <- function(d, bn, sn){
        mean(d[bn, grep(sn, colnames(d))], na.rm = TRUE)
      }
      br <- sapply(names(reference), thisMean, d = d, sn = spec[1])
      br <- names(br)[which.min(br)]
    }
    reference <- reference[br, ]
    br <- paste("REF", br, sep = "_")
    rownames(reference) <- br
  } else {
    br  <- "^REF_"
  }
  
  # profile alignment is far quicker
  obj <- mafft.merge(list(reference, obj), 
                     exec = megProj@align.exe,
                     thread = 1)
  
  ## handle reverse complements
  ## --------------------------
  rc <- grep("^_R_", rownames(obj))
  if (length(rc)){
    spec.id <- grep(spec[1], rownames(obj))
    rownames(obj) <- gsub("^_R_", "", rownames(obj))
    ## realignment if only subset of sequences concerned
    if ( length(spec.id) > length(rc) ){
      gt <- splitGiTaxon(rownames(obj)[rc], sep = " ")
      rc <- paste("UPDATE", acc.tab,
                  "SET status = 'aligned-RC'",
                  "WHERE", sql.wrap(gt[, 1], term = "gi", BOOL = NULL),
                  "AND", sql.wrap(gt[, 2], term = "taxon", BOOL = NULL))
      dbSendQuery(conn, rc)
      rc <- mafft(obj[spec.id, ], exec = megProj@align.exe)
      dbWriteDNA(conn, acc.tab, rc)
      
      ## all sequences concerned
    } else {
      rc <- obj[rc, ]
      rc <- deleteEmptyCells(rc, TRUE)
      dbWriteDNA(conn, acc.tab, rc, status = "aligned-RC")
    } 
  }
  
  ## crop sequences + update database
  ## --------------------------------
  n <- ncol(obj)
  obj <- cropToReference(obj)
  if (ncol(obj) < n){
    
    ## identify sequences that do not overlap with reference
    ## and mark them as 'excluded (outside reference)'
    id <- apply(obj, 1, 
                function(x) ifelse(all(unique(x) %in% as.raw(c(240, 2, 4))), 
                                   FALSE, TRUE))
    id <- names(which(!id))
    if ( length(id) > 0 ){
      id <- data.frame(splitGiTaxon(id), sep = " ")
      SQL <- paste("UPDATE", acc.tab, 
                   "SET status = 'excluded (outside reference)'",
                   "WHERE", wrapSQL(id$gi, "gi", "=", NULL), 
                   "AND", wrapSQL(id$taxon, "taxon", "=", NULL))
      lapply(SQL, dbSendQuery, conn = conn)
      ## delete non-overlapping sequences from alignment
      obj <- deleteEmptyCells(obj, quiet = TRUE)
    }
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
  d <- data.frame(splitGiTaxon(names(d), sep = " "), 
                  identity = d,
                  coverage = cv)
  SQL <- paste("UPDATE", acc.tab, 
               "SET", wrapSQL(d$identity, "identity", "=", NULL), 
               ",", wrapSQL(d$coverage, term = "coverage", "=", NULL),
               "WHERE", wrapSQL(d$gi, term = "gi", "=", NULL), 
               "AND", wrapSQL(d$taxon, "taxon", "=", NULL))
  lapply(SQL, dbSendQuery, conn = conn)
  dbDisconnect(conn)
} 