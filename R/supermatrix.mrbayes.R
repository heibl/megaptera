## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-02-17)

supermatrix.mrbayes <- function(megapteraProj, min.n.seq, masked = TRUE, 
                        partition,
                        coverage.locus = 0.5,
                        exclude.locus, exclude.species,
                        core.locus, core.species,
                        squeeze.outgroup, 
                        lset, mcmc){
  
  ## join taxonomy and locus tables
  ## ------------------------------
  conn <- dbconnect(megapteraProj)
  tab <- "SELECT * FROM taxonomy INNER JOIN locus USING (spec)"
  tab <- dbGetQuery(conn, tab)
  tab[is.na(tab)] <- "excluded"
  cols <- grep("_blocks", names(tab))
  tab <- tab[apply(tab[, cols], 1, 
                   function(x) length(grep("excluded", x)) != length(x)), ]
  
  ## outgroup
  ## --------
  outgroup <- megapteraProj@taxon@outgroup
  if ( is.list(outgroup) ) outgroup <- sapply(outgroup, head, 1)
  if ( !unique(is.Linnean(outgroup)) )
    outgroup <- tab$spec[which(tab == outgroup, arr.ind = TRUE)[, "row"]]
  outgroup <- gsub(" ", "_", outgroup)
  
  ## determine tables that contain more than min.n.seq species
  ## ---------------------------------------------------------
  if ( missing(min.n.seq) ){
    tabnames <- names(tab[cols])
  } else {
    tabnames <- apply(tab[, cols], 2, function(x) length(grep("selected", x)))
    tabnames <- names(tabnames[tabnames >= min.n.seq])
  }
  tabnames <- paste("spec", gsub("_blocks", "", tabnames), sep = "_")
  tabnames <- gsub("__", "_", tabnames)
  
  ## exclude loci from concatenation
  ## -------------------------------
  if ( !missing(exclude.locus) ){
    cat("\n.. excluding locus:", exclude.locus, " ..")
    exclude.locus <- paste(exclude.locus, collapse = "|")
    tabnames <- tabnames[-grep(exclude.locus, tabnames)]
  }
  
  ## select and read individual alignments
  ## -------------------------------------
  cat("\n.. reading", length(tabnames), "alignments ..")
  x <- lapply(tabnames, dbReadDNA, conn = conn, masked = masked)
  test <- sapply(x, function(a) ncol(a) == 2)
  if ( any(test) ){
    stop("no sequences found in '", tabnames[test], "' - have you run stepG and stepH?")
  }
  x <- lapply(x, trimEnds)
  dbDisconnect(conn)
  names(x) <- gsub("^spec_", "", tabnames)
  spec.set <- lapply(x, rownames)
  spec.set <- table(unlist(spec.set))
  spec.set <- data.frame(spec = names(spec.set),
                         freq = spec.set,
                         stringsAsFactors = FALSE)
  nspec <- nrow(spec.set)
  
  ## exclude species from alignments that have
  ## less than 'coverage.locus' percent sites
  ## ----------------------------------------
  cat("\n.. excluding snippets (<", coverage.locus, "% sites) ..")
  exclude.snippets <- function(a, coverage.locus){
    cv <- coverage(a)
    cv <- names(cv)[cv >= coverage.locus]
    a[cv, ]
  }
  x <- lapply(x, exclude.snippets, coverage.locus = coverage.locus)
  ## any species lost?
  spec.set2 <- lapply(x, rownames)
  spec.set2 <- table(unlist(spec.set2))
  spec.set2 <- data.frame(spec = names(spec.set2),
                          freq = spec.set2,
                          stringsAsFactors = FALSE)
  lost <- setdiff(spec.set$spec, spec.set2$spec)
  if ( length(lost) > 0 ){
    cat("\n   WARNING:", length(lost), "species lost:",
        paste("\n  -", lost))
  }
  
  ## core set of species
  ## -------------------
  if ( !missing(core.locus) ){
    cat("\n.. making core dataset ..")
    core.species <- lapply(x[core.locus], function(x) rownames(x))
    core.species <- unique(unlist(core.species))
    subset.alignment <- function(a, s) deleteEmptyCells(a[rownames(a) %in% s, ], quiet = TRUE)
    x <- lapply(x, subset.alignment, s = core.species)
  }
  
  ## create partitions
  ## -----------------
  if ( !missing(partition) ){
    partition <- lapply(partition, intersect, names(x))
    x <- x[match(unlist(partition), names(x))]
    concatenate <- function(dna, loci){
      dna <- dna[loci]
      do.call(cbind.DNAbin, c(dna, fill.with.gaps = TRUE))
    }
    x <- xx <- lapply(partition, concatenate, dna = x)
  } else {
    xx <- x
  }
  ngene <- length(xx)
  p <- cbind(rep(1, length(x)), sapply(x, ncol))
  for ( i in 2:nrow(p) ){
    p[i, 1] <- p[i - 1, 2] + 1
    p[i, 2] <- p[i, 1] + p[i, 2] -1
  }
  ## partitions in RAxML format:
  p <- paste("DNA, ", rownames(p), " = ", p[, 1], "-", p[, 2], sep = "")
  
  ## create SUPERMATRIX
  ## ------------------
  cat("\n.. concatenating alignments ..")
  x <- do.call(cbind.DNAbin, c(x, fill.with.gaps = TRUE))
  
  ## exclude species by user decision
  ## --------------------------------
  if ( !missing(exclude.species) ){
    exclude.species <- intersect(exclude.species, rownames(x))
    nes <- length(exclude.species)
    if ( nes > 0 ){
      cat("\n.. excluding ", nes ," (", 
          round(nes/nrow(x), 2), 
          "%) species by user decision ..", sep = "")
      x <- x[!rownames(x) %in% exclude.species, ]
    }
  }
  
  ## outgroup
  ## ---------
  outgroup <- intersect(outgroup, rownames(x))
  if ( !missing(squeeze.outgroup) ){
    cat("\n.. creating denser outgroup ..")
    o <- deleteEmptyCells(x[outgroup, ], quiet = TRUE)
    nn <- as.raw(c(240, 2, 4))
    names(nn) <- c("n", "?", "-")
    nn <- apply(o, 1, function(x, n) length(which((x %in% nn))), n = nn)
    nn <- sort(nn, decreasing = TRUE)
    outgroup <- names(tail(nn, squeeze.outgroup))
    nn <- head(nn, -squeeze.outgroup)
    x <- x[!rownames(x) %in% names(nn), ]
    x <- deleteEmptyCells(x, quiet = TRUE)
  }
  
  ingroup <- setdiff(rownames(x), outgroup)
  
  ## write NEXUS file for MrBayes
  ## ----------------------------
  cat("\n.. writing NEXUS file for MrBayes ..")
  masked <- ifelse(masked, "masked", "")
  fn <- paste("supermatrix", nrow(x), ngene, ncol(x), sep = "-")
  fn <- paste(fn, "bayes.nex", sep = ".")
  # rownames(x) <- gsub("-", "_", rownames(x)) comment tentative!
  
  prset <- mrbayes.prset(topologypr = "constraint(ingroup)",
                         brlenspr = "clock:uniform",
                         clockvarpr = "igr",
                         igrvarpr = "exponential(10.00)")
  unlink <- list(statefreq = "all", revmat = "all",
                 shape = "all", pinvar = "all")
  obj <- mrbayes(xx, fn, constraint = list(ingroup = ingroup),
                 unlink = unlink,
                 prset = prset, lset = lset, mcmc = mcmc)
}

