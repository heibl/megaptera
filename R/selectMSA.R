## This code is part of the megaptera package
## © C. Heibl 2014 (last update 20120-02-27)

#' @title Selection of Multiple Sequence Alignments
#' @importFrom ape cbind.DNAbin
#' @export

selectMSA <- function(megProj, min.n.seq = 3, blocks = "split", 
                      row.confid = 0, col.confid = 0,
                      trim.ends = 0,
                      locus.coverage = 0.0,
                      global.coverage = 0.0,
                      subset.locus, subset.species,
                      exclude.locus, exclude.species,
                      core.locus, core.species,
                      best.sampled.congeneric = FALSE,
                      equal.taxon.sets = TRUE,
                      use.outgroup = TRUE,
                      protect.outgroup = FALSE, squeeze.outgroup){
  
  ## INITIAL CHECKS + ADJUSTMENTS
  ## ----------------------------
  if (!inherits(megProj, "megapteraProj")) stop("'megProj' is not of class 'megapteraProj'")
  tip.rank <- megProj@taxon@tip.rank
  if (length(min.n.seq) == 1) min.n.seq <- rep(min.n.seq, 2)
  blocks <- match.arg(blocks, c("concatenate", "ignore", "split"))
  if (blocks == "concatenate") blocks <- "split"
  
  ## Get outgroup
  ## ------------
  outgroup <- dbReadTaxonomy(megProj)
  outgroup <- lapply(megProj@taxon@outgroup, taxdumpChildren,
                     tax = outgroup, tip.rank = "species")
  outgroup <- gsub(" ", "_", do.call(rbind, outgroup)$taxon)
  
  #############################################
  ##  PART A: Determine tables (loci) to import
  #############################################
  
  ## A1: Get list of available tables
  ## --------------------------------
  cat("Looking for database tables ... ")
  tabs <- checkBlocks(megProj, plot = FALSE, subset = subset.species)
  # tabs <- checkBlocks(megProj, plot = FALSE)
  cat(length(tabs), " found:", paste("\n-", sort(names(tabs))), sep = "")
  
  ## A2: Determine tables that contain more than min.n.seq species
  ## -------------------------------------------------------------
  id <- sapply(tabs, function(z) any(z >= min.n.seq[1]))
  tabs <- names(tabs)[id]
  cat("\n", length(tabs), " of these contain at least ", min.n.seq[1], 
      " species:", sort(paste("\n-", tabs)), sep = "")
  
  ## A3: User-defined subset loci (optional)
  ## ---------------------------------------
  if (!missing(subset.locus)){
    cat("\nSubsetting to loci:", subset.locus)
    subset.locus <- paste(subset.locus, collapse = "|")
    tabs <- tabs[grep(subset.locus, tabs)]
  }
  
  ## A4: User-defined exclusion of loci (optional)
  ## ---------------------------------------------
  if (!missing(exclude.locus)){
    cat("\nExcluding locus:", exclude.locus)
    exclude.locus <- paste(exclude.locus, collapse = "|")
    tabs <- tabs[-grep(exclude.locus, tabs)]
  }
  
  #############################################
  ##  PART B: Import alignments
  #############################################
  
  ## B1: Select and read individual alignments
  ## -----------------------------------------
  cat("\nReading", length(tabs), "alignments ... ")
  if (missing(subset.species)){
    x <- lapply(tabs, dbReadMSA, x = megProj, blocks = blocks,
                label = "taxon", confid.scores = "col.means", 
                row.confid = row.confid, col.confid = col.confid)
  } else {
    x <- lapply(tabs, dbReadMSA, x = megProj, taxon = subset.species,
                label = "taxon",
                regex = FALSE, blocks = blocks, 
                row.confid = row.confid, col.confid = col.confid)
  }
  names(x) <- tabs
  cat("OK")
  if (any(sapply(x, is.null))) x <- x[!sapply(x, is.null)]
  
  ## B2: Exclude species
  ## -------------------
  if (!missing(exclude.species)){
    cat("\nExcluding species locus-wise (user-decision)")
    for (loci in names(exclude.species)){
      a <- x[[loci]]
      cs <- attr(a, "cs")
      id <- rownames(a) %in% gsub(" ", "_", exclude.species[[loci]])
      cat("\n- ", loci, ": ", paste(sort(exclude.species[[loci]]), collapse = (", ")), sep = "")
      a <- a[!id, ]
      attr(a, "cs") <- cs
      x[[loci]] <- a
    }
  }
  
  ## B3: Handle blocks
  ## -----------------
  block.id <- which(sapply(x, is.list))
  cat("\n", length(block.id), " alignments are split into blocks", sep = "")
  if (length(block.id)){
    concatenateBlocks <- function(ali, n){
      ali <- ali[which(sapply(ali, nrow) >= n)]
      if (length(ali)){
        ali <- do.call(cbind.DNAbin, c(ali, fill.with.gaps = TRUE))
      } else {
        ali <- ali[[1]]
      }
      ali
    }
    cat("\nKeeping only blocks >", min.n.seq[2], "sequence")
    x[block.id] <- lapply(x[block.id], concatenateBlocks, n = min.n.seq[2])
  }
  
  ## B4: Trim tapering ends of alignment
  ## -----------------------------------
  cat("\nTrimming alignment ends (columns with <", 
      round(trim.ends * 100, 1), "% sequence information) ... ")
  n <- sapply(x, ncol)
  x <- lapply(x, trimEnds, min.n.seq = trim.ends)
  cat(sum(n - sapply(x, ncol)), "columns deleted")
  
  names(x) <- gsub(paste0("^", tip.rank, "_"), "", tabs)
  spec.set <- lapply(x, rownames)
  spec.set <- table(unlist(spec.set))
  spec.set <- data.frame(spec = names(spec.set),
                         freq = spec.set,
                         stringsAsFactors = FALSE)
  nspec <- nrow(spec.set)
  
  ## B5: Exclude species from alignments that have
  ## less than 'locus.coverage' percent sites
  ## ----------------------------------------
  if (locus.coverage > 0){
    cat("\nExcluding snippets (<", locus.coverage, "% sites) ... ")
    exclude.snippets <- function(a, locus.coverage){
      id <- coverage(a) >= locus.coverage
      cs <- attr(a, "cs")
      a <- a[id, ]
      if (!is.null(cs)){
        if (is.matrix(cs)){
          cs <- cs[id,]
        } else {
          if (!all(id)) warning("change in alignment not traceable in confidence scores")
        }
        attr(a, "cs") <- cs
      }
      a
    }
    x <- lapply(x, exclude.snippets, locus.coverage = locus.coverage)
  }
  
  ## B5: Core set of species (optional)
  ## ----------------------------------
  if (!missing(core.locus)){
    cat("\nMaking core dataset ... ")
    core.species <- lapply(x[core.locus], function(x) rownames(x))
    core.species <- unique(unlist(core.species))
    if (protect.outgroup) core.species <- union(core.species, outgroup)
    subset.alignment <- function(a, s) deleteEmptyCells(a[rownames(a) %in% s, ], quiet = TRUE)
    x <- lapply(x, subset.alignment, s = core.species)
    x <- lapply(x, deleteSpecies, species = nn)
    cat("OK")
    
    ## Any species lost?
    spec.set2 <- lapply(x, rownames)
    spec.set2 <- table(unlist(spec.set2))
    spec.set2 <- data.frame(spec = names(spec.set2),
                            freq = spec.set2,
                            stringsAsFactors = FALSE)
    lost <- setdiff(spec.set$spec, spec.set2$spec)
    nb.lost <- length(lost)
    if (nb.lost){
      if (nb.lost > 12) lost <- c(head(lost), paste0("[", nb.lost - 12, " sequences]"), tail(lost))
      cat("WARNING:", nb.lost, "species lost:",
          paste("\n-", lost))
    } else {
      cat("OK")
    }
  }
  
  ## Delete outgroup
  ## ---------------
  if (!use.outgroup){
    cat("\nPreparing alignments without outgroup ... ")
    x <- lapply(x, deleteSpecies, delete = outgroup)
    cat("OK")
  }
  
  ## Create denser outgroup
  ## ----------------------
  if (!missing(squeeze.outgroup)){
    cat("\nCreating denser outgroup ... ")
    o <- do.call(cbind.DNAbin, c(x, fill.with.gaps = TRUE))
    o <- deleteEmptyCells(o[outgroup, ], quiet = TRUE)
    nn <- as.raw(c(240, 2, 4))
    names(nn) <- c("n", "?", "-")
    nn <- apply(o, 1, function(x, n) length(which((x %in% nn))), n = nn)
    nn <- sort(nn, decreasing = TRUE)
    outgroup <- names(tail(nn, squeeze.outgroup))
    x <- lapply(x, deleteSpecies, delete = head(nn, -squeeze.outgroup))
    cat("OK")
  }
  
  ## Create partition matrix
  ## -----------------------
  # if (!missing(partition)){
  #   partition <- lapply(partition, intersect, names(x))
  #   x <- x[match(unlist(partition), names(x))]
  #   concatenate <- function(dna, loci){
  #     dna <- dna[loci]
  #     do.call(cbind.DNAbin, c(dna, fill.with.gaps = TRUE))
  #   }
  #   x <- xx <- lapply(partition, concatenate, dna = x)
  # } else {
  #   xx <- x
  # }
  p <- cbind(rep(1, length(x)), sapply(x, ncol))
  if (nrow(p) > 1){
    for (i in 2:nrow(p)){
      p[i, 1] <- p[i - 1, 2] + 1
      p[i, 2] <- p[i, 1] + p[i, 2] -1
    }
  }
  ## partitions in RAxML format:
  # p <- paste0("DNA, ", rownames(p), " = ", p[, 1], "-", p[, 2])
  
  ## Use cbind.DNAbin to fill missing data and then split again if
  ## equal taxon sets are requires (e.g. for linking trees in BEAST)
  ## ---------------------------------------------------------------
  if (equal.taxon.sets){
    x <- do.call(cbind.DNAbin, c(x, fill.with.gaps = TRUE))
    x <- apply(p, 1, function(x, z) x[, z[1]:z[2]], x = x)
  }
  
  x
}

## Function to delete species from MSA with CS
## -------------------------------------------
deleteSpecies <- function(msa, delete, keep){
  
  keep <- setdiff(rownames(msa), delete)
  cs <- attr(msa, "cs")
  ## if present cs is considered to be a list
  if (is.matrix(cs)) stop("debug me! [1]")
  msa <- msa[rownames(msa) %in% keep, ]
  id <- identifyEmptyCells(msa, quiet = TRUE)
  if (!is.null(cs)){
    if (length(id$col)){
      msa <- msa[, -(id$col)]
      cs <- cs[-(id$col)]
    } 
    if (ncol(msa) != length(cs)) stop("debug me! [2]")
    attr(msa, "cs") <- cs
  }
  msa
}