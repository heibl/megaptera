## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2021-03-15)

#' @title Selection of Multiple Sequence Alignments
#' @importFrom ape cbind.DNAbin
#' @importFrom crayon %+% bold cyan magenta red silver
#' @export

selectMSA <- function(megProj, min.n.seq = 3, blocks = "split", 
                      row.confid = 0, col.confid = 0,
                      trim.ends = 0,
                      locus.coverage = 0.0,
                      global.coverage = 0.0,
                      subset.locus = NULL, subset.species = NULL,
                      exclude.locus = NULL, exclude.species = NULL,
                      core.locus = NULL, core.species = NULL,
                      best.sampled.congeneric = FALSE,
                      equal.taxon.sets = TRUE,
                      use.outgroup = TRUE, protect.outgroup = FALSE, squeeze.outgroup = NULL){
  
  ## Default argument values for debugging
  ## min.n.seq = 3; blocks = "split"; row.confid = 0; col.confid = 0
  ## subset.locus <- subset.species <- exclude.locus <- exclude.species <- core.locus <- core.species <- NULL
  
  ## INITIAL CHECKS + ADJUSTMENTS
  ## ----------------------------
  if (!inherits(megProj, "megapteraProj")) stop("'megProj' is not of class 'megapteraProj'")
  tip.rank <- megProj@taxon@tip.rank
  if (length(min.n.seq) == 1) min.n.seq <- rep(min.n.seq, 2)
  blocks <- match.arg(blocks, c("concatenate", "ignore", "split"))
  if (blocks == "concatenate") blocks <- "split"
  
  ## Get ingroup and outgroup
  ## ------------------------
  tax <- dbReadTaxonomy(megProj)
  ingroup <- lapply(megProj@taxon@ingroup, taxdumpChildren,
                     tax = tax, tip.rank = "species")
  ingroup <- gsub(" ", "_", do.call(rbind, ingroup)$taxon)
  og <- outgroup(megProj, sep = "_")
  
  #############################################
  ##  PART A: Determine tables (loci) to import
  #############################################
  
  ## A1: Get list of available tables
  ## --------------------------------
  cat(silver("Checking available loci ... "))
  tabs <- checkBlocks(megProj, plot = FALSE, subset = subset.species)
  cat(length(tabs), " found:\n", paste(" > ", sort(names(tabs)), "\n"), sep = "")
  
  ## A2: Determine tables that contain more than min.n.seq species
  ## -------------------------------------------------------------
  id <- sapply(tabs, function(z) any(z >= min.n.seq[1]))
  tabs <- names(tabs)[id]
  cat(silver(length(tabs), " of these contain at least ", min.n.seq[1], 
      " species:", sort(paste(" > ", tabs, "\n"))), sep = "")
  
  ## A3: User-defined subset loci (optional)
  ## ---------------------------------------
  if (!is.null(subset.locus)){
    cat("\nSubsetting to loci:", subset.locus)
    subset.locus <- paste(subset.locus, collapse = "|")
    tabs <- tabs[grep(subset.locus, tabs)]
  }
  
  ## A4: User-defined exclusion of loci (optional)
  ## ---------------------------------------------
  if (!is.null(exclude.locus)){
    cat("\nExcluding locus:", exclude.locus)
    exclude.locus <- paste(exclude.locus, collapse = "|")
    tabs <- tabs[-grep(exclude.locus, tabs)]
  }
  
  #############################################
  ##  PART B: Import alignments
  #############################################
  
  ## B1: Select and read individual alignments
  ## -----------------------------------------
  cat(silver("Reading " %+% magenta$bold(length(tabs)) %+% " alignments ... "))
  if (is.null(subset.species)){
    obj <- lapply(tabs, dbReadMSA, x = megProj, blocks = blocks,
                label = "taxon", confid.scores = "col.means", 
                row.confid = row.confid, col.confid = col.confid)
  } else {
    obj <- lapply(tabs, dbReadMSA, x = megProj, taxon = subset.species,
                label = "taxon",
                regex = FALSE, blocks = blocks, 
                row.confid = row.confid, col.confid = col.confid)
  }
  names(obj) <- tabs
  if (any(sapply(obj, is.null))) obj <- obj[!sapply(obj, is.null)]
  cat(green("OK\n"))
  
  ## taxonomic coverage
  taxCov(obj, ingroup, og)
 
  
  
  ## B2: Exclude species
  ## -------------------
  if (!is.null(exclude.species)){
    cat("Excluding species locus-wise (user-decision)")
    for (loci in names(exclude.species)){
      a <- obj[[loci]]
      cs <- attr(a, "cs")
      id <- rownames(a) %in% gsub(" ", "_", exclude.species[[loci]])
      cat("\n- ", loci, ": ", paste(sort(exclude.species[[loci]]), collapse = (", ")), sep = "")
      a <- a[!id, ]
      attr(a, "cs") <- cs
      obj[[loci]] <- a
    }
  }
  
  ## B3: Handle blocks
  ## -----------------
  block.id <- which(sapply(obj, is.list))
  cat(silver(length(block.id), " alignments are split into blocks\n"), sep = "")
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
    cat("Keeping only blocks >", min.n.seq[2], "sequence\n")
    obj[block.id] <- lapply(obj[block.id], concatenateBlocks, n = min.n.seq[2])
  }
  
  ## B4: Trim tapering ends of alignment
  ## -----------------------------------
  cat("Trimming alignment ends (columns with <", 
      round(trim.ends * 100, 2), "% sequence information) ... ")
  n <- sapply(obj, ncol)
  obj <- lapply(obj, trimEnds, min.n.seq = trim.ends)
  cat(sum(n - sapply(obj, ncol)), "columns deleted")
  
  names(obj) <- gsub(paste0("^", tip.rank, "_"), "", tabs)
  spec.set <- lapply(obj, rownames)
  spec.set <- table(unlist(spec.set))
  spec.set <- data.frame(spec = names(spec.set),
                         freq = spec.set,
                         stringsAsFactors = FALSE)
  nspec <- nrow(spec.set)
  cat(green("OK\n"))
  ## taxonomic coverage
  taxCov(obj, ingroup, og)
  
  ## B5: Exclude species from alignments that have
  ## less than 'locus.coverage' percent sites
  ## ----------------------------------------
  if (locus.coverage > 0){
    cat(silver("Excluding snippets (<", locus.coverage, "% sites) ... "))
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
    obj <- lapply(obj, exclude.snippets, locus.coverage = locus.coverage)
    cat(green("OK\n"))
    ## taxonomic coverage
    taxCov(obj, ingroup, og)
  }
  
  ## B5: Core set of species (optional)
  ## ----------------------------------
  if (!is.null(core.locus)){
    cat(silver("Making core dataset ... "))
    core.species <- lapply(obj[core.locus], function(x) rownames(x))
    core.species <- unique(unlist(core.species))
    if (protect.outgroup) core.species <- union(core.species, outgroup)
    subset.alignment <- function(a, s) deleteEmptyCells(a[rownames(a) %in% s, ], quiet = TRUE)
    obj <- lapply(obj, subset.alignment, s = core.species)
    obj <- lapply(obj, deleteSpecies, species = nn)
    cat(green("OK\n"))
    
    ## Any species lost?
    spec.set2 <- lapply(obj, rownames)
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
      cat(green("OK\n"))
    }
    ## taxonomic coverage
    taxCov(obj, ingroup, og)
  }
  
  ## Delete outgroup
  ## ---------------
  if (!use.outgroup){
    cat("\nPreparing alignments without outgroup ... ")
    obj <- lapply(obj, deleteSpecies, delete = og)
    cat("OK")
  }
  
  ## Create denser outgroup
  ## ----------------------
  if (!is.null(squeeze.outgroup)){
    cat("\nCreating denser outgroup ... ")
    o <- do.call(cbind.DNAbin, c(obj, fill.with.gaps = TRUE))
    o <- deleteEmptyCells(o[og, ], quiet = TRUE)
    nn <- as.raw(c(240, 2, 4))
    names(nn) <- c("n", "?", "-")
    nn <- apply(o, 1, function(z, n) length(which((z %in% nn))), n = nn)
    nn <- sort(nn, decreasing = TRUE)
    og <- names(tail(nn, squeeze.outgroup))
    obj <- lapply(obj, deleteSpecies, delete = head(nn, -squeeze.outgroup))
    cat("OK")
  }
  
  ## Create partition matrix
  ## -----------------------
  # if (!missing(partition)){
  #   partition <- lapply(partition, intersect, names(obj))
  #   obj <- obj[match(unlist(partition), names(obj))]
  #   concatenate <- function(dna, loci){
  #     dna <- dna[loci]
  #     do.call(cbind.DNAbin, c(dna, fill.with.gaps = TRUE))
  #   }
  #   obj <- xx <- lapply(partition, concatenate, dna = obj)
  # } else {
  #   xx <- obj
  # }
  p <- cbind(rep(1, length(obj)), sapply(obj, ncol))
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
    stop("debug me!")
    obj <- do.call(cbind.DNAbin, c(obj, fill.with.gaps = TRUE))
    obj <- apply(p, 1, function(x, z) x[, z[1]:z[2]], x = obj)
  }
  
  obj
}

## Function to delete species from MSA with CS
## -------------------------------------------
## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2021-03-12)

#' @importFrom ips identifyEmptyCells
#' @export

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