## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-04-08)

# TO DO: marker-wise deletion of species (line 38)

supermatrix <- function(megapteraProj, min.n.seq = 3, 
                        masked = TRUE, 
                        partition,
                        subset.locus, coverage.locus = 0.5,
                        exclude.locus, exclude.species,
                        core.locus, core.species,
                        squeeze.outgroup){
  
  tip.rank <- megapteraProj@taxon@tip.rank
  
  ## determine tables that contain more than min.n.seq species
  ## ---------------------------------------------------------
  tabs <- checkBlocks(megapteraProj, plot = FALSE)
  id <- sapply(tabs, function(z) any(z >= min.n.seq))
  tabs <- names(tabs)[id]
  
  ## subset loci to locus
  ## ---------------------
  if ( !missing(subset.locus) ){
    cat("\n.. subsetting to loci:", subset.locus, " ..")
    subset.locus <- paste(subset.locus, collapse = "|")
    tabs <- tabs[grep(subset.locus, tabs)]
  }
  
  ## exclude loci from concatenation
  ## -------------------------------
  if ( !missing(exclude.locus) ){
    cat("\n.. excluding locus:", exclude.locus, " ..")
    exclude.locus <- paste(exclude.locus, collapse = "|")
    tabs <- tabs[-grep(exclude.locus, tabs)]
  }
  
  ## select and read individual alignments
  ## -------------------------------------
  cat("\n.. reading", length(tabs), "alignments ..")
  conn <- dbconnect(megapteraProj)
  x <- lapply(tabs, dbReadDNA, 
              x = megapteraProj, taxon = ".*", regex = TRUE,
              ignore.excluded = TRUE,
              blocks = "split",
              masked = masked)
#   test <- sapply(x, function(a) ncol(a) == 2)
#   if ( any(test) ){
#     stop("no sequences found in '", tabnames[test], "' - have you run stepG and stepH?")
#   }
  ## handle blocks
  ## -------------
  block.id <- which(sapply(x, is.list))
  cat("\n..", length(block.id), "alignments are split into blocks ..")
  concatenateBlocks <- function(ali, n){
    ali <- ali[which(sapply(ali, nrow) >= n)]
    if ( length(ali) > 1 ){
      ali <- do.call(cbind.DNAbin, c(ali, fill.with.gaps = TRUE))
    } else {
      ali <- ali[[1]]
    }
    ali
  }
  x[block.id] <- lapply(x[block.id], concatenateBlocks, n = min.n.seq)
  
  x <- lapply(x, trimEnds, min.n.seq = min.n.seq)
  
  names(x) <- gsub(paste("^", tip.rank, "_", sep = ""), "", tabs)
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
  ## --------
  outgroup <- paste("SELECT", tip.rank, 
                    "FROM taxonomy WHERE tag ~ 'outgroup'")
  outgroup <- dbGetQuery(conn, outgroup)[, tip.rank]
  outgroup <- intersect(outgroup, rownames(x))
  cat("\n.. number of available outgroup species:", length(outgroup),  "..")
  dbDisconnect(conn)
  
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
  
  ## make filenames (from here on 'x' does not change any more)
  ## ----------------------------------------------------------
  masked <- ifelse(masked, "masked", "")
  fn <- paste("supermatrix", nrow(x), ngene, ncol(x), sep = "-")
  ext <- c("tre", "phy", "nex", "partitions", "outgroup")
  fns <- paste(fn, ext, sep = ".")
  names(fns) <- ext
  
  
  ## prepare guide tree 
  ## ------------------
  cat("\n.. preparing comprehensive guidetree ..")
  if ( inherits(megapteraProj@taxon, "taxonGuidetree") ){
    gt <- comprehensiveGuidetree(megapteraProj, 
                                 tip.rank = tip.rank,
                                 subset = x)
  } else { 
    tax <- dbReadTaxonomy(megapteraProj, subset = x)
    gt <- tax2tree(tax, tip.rank = tip.rank)
  }
  gt <- drop.tip(gt, setdiff(gt$tip.label, rownames(x)))
  
  ## write outgroup + partitions files
  ## ---------------------------------
  og <- paste(outgroup, collapse = ",")
  clip <- pipe("pbcopy", "w")
  write(og, file = clip)
  close(clip)
  write(og, fns["outgroup"])
  write(p, fns["partitions"])
  
  ## write data as PHY and NEX
  ## -------------------------
  cat("\n.. writing supermatrix to files ..")
  write.tree(gt, fns["tre"])
  write.phy(x, fns["phy"])
  rownames(x) <- gsub("-", "_", rownames(x))
  write.nex(x, fns["nex"])
  
  ## zip and return
  ## --------------
  cat("\n.. zipping supermatrix ..\n")
  zip(zipfile = fn, files = fns[c("phy", "partitions", "outgroup")])
  obj <- list(zipname = fn,
              supermatrix = x,
              outgroup = og,
              partitions = p,
              guide.tree = gt)
  obj
}

