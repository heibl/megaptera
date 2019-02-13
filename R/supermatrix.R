## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2019-02-06)

# TO DO: marker-wise deletion of species (line 38)

#' @title Retrieval and Concatentation of Multiple Sequence Alignments
#' @description Concatenates multiple sequence alignments of individual loci
#'   into a supermatrix.
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @param min.n.seq Numeric, the minimum number of sequences in any alignment
#'   that is required to be included in the supermatrix.
#' @param reliability Numeric between 0 and 1, giving the minimum reliability
#'   score for a column in an alignment.
#' @param blocks A character string indicating how to handle alignment blocks:
#'   \code{"split"} causes blocks to be returned as elements of a list.
#'   \code{"concatenate"} means, blocks will be concatendated and returned as a
#'   single alignment.
#' @param partition A named list, each element of which contains the names of
#'   the loci to group together in one partition. The default is to assign each
#'   locus to its own partition. Partitioning of first/second and third
#'   nucleotides positions is not yet possible, but might be implemented in the
#'   future.
#' @param trim.ends A \code{numeric} giving the required minimum number of
#'   sequences having an non-ambiguous base character (a, c, g, t) in the first
#'   and last position of the alignment; defaults to \code{0}, which means no
#'   trimming. topology. Can also be given as a fraction.
#' @param coverage.locus Numeric between 0 and 1 giving the required minimum
#'   coverage of any species in any alignment.
#' @param subset.locus A vector of mode \code{"character"} for choosing a subset
#'   of loci from the loci available.
#' @param subset.species A vector of mode \code{"character"} for choosing a
#'   subset of the species from the total species available.
#' @param exclude.locus A vector of mode \code{"character"} giving the names of
#'   the loci in the database that will be excluded from the concatenation.
#' @param exclude.species \emph{Currently unused.}
#' @param core.locus A vector of mode \code{"character"} giving the names of the
#'   'core' loci: The resulting supermatrix will only contain species that are
#'   contained in the 'core' loci. This option is intended to create denser
#'   supermatrices.
#' @param core.species \emph{Currently unused.}
#' @param best.sampled.congeneric Logical, keep all but the best-sampled species
#'   in every genus.
#' @param squeeze.outgroup Numeric, can be given to reduce the number of
#'   outgroup species: The function will select the \code{squeeze.outgroup}
#'   outgroup species with the best coverage. The idea is to create a more
#'   densely sampled outgroup.
#' @return a list with four elements: \item{zipname}{a character string giving
#'   the name of the files produced.} \item{supermatrix}{a matrix of class
#'   \code{\link{DNAbin}}.} \item{outgroup}{a vector of mode \code{"character"}
#'   giving the names of the species used as an outgroup.} \item{partitions}{a
#'   vector of mode \code{"character"} giving the partitions of the supermatrix
#'   in the same spelling as accepted by RAxML.}
#'
#'   In addition, three files are written to the working directory: (1) a
#'   NEXUS-formatted file and (2) a PHYLIP-formatted file containing the
#'   sequence matrix and (3) a zipped directory containing the PHYLIP-formatted
#'   sequence matrix plus output and partitions as separate ASCII-formatted
#'   files.
#' @importFrom ape write.tree
#' @importFrom utils zip
#' @export

supermatrix <- function(megProj, min.n.seq = 3, 
                        reliability = 0, blocks = "split", 
                        partition, trim.ends = 0,
                        coverage.locus = 0.5,
                        subset.locus, subset.species,
                        exclude.locus, exclude.species,
                        core.locus, core.species,
                        best.sampled.congeneric = FALSE,
                        squeeze.outgroup){
  
  ## INITIAL CHECKS + ADJUSTMENTS
  ## ----------------------------
  if (!inherits(megProj, "megapteraProj")) stop("'megProj' is not of class 'megapteraProj'")
  tip.rank <- megProj@taxon@tip.rank
  if (length(min.n.seq) == 1) min.n.seq <- rep(min.n.seq, 2)
  blocks <- match.arg(blocks, c("concatenate", "ignore", "split"))
  if (blocks == "concatenate") blocks <- "split"
  
  #############################################
  ##  PART 1: Determine tables (loci) to import
  #############################################
  
  ## Determine tables that contain more than min.n.seq species
  ## ---------------------------------------------------------
  cat("Looking for database tables ... ")
  tabs <- checkBlocks(megProj, plot = FALSE, subset = subset.species)
  cat(length(tabs), " found:", paste("\n-", names(tabs)), sep = "")
  id <- sapply(tabs, function(z) any(z >= min.n.seq[1]))
  tabs <- names(tabs)[id]
  cat("\n", length(tabs), " of these contain at least ", min.n.seq[1], " species:", 
      paste("\n-", tabs), sep = "")
  
  ## User-defined subset loci (optional)
  ## -----------------------------------
  if (!missing(subset.locus)){
    cat("\nSubsetting to loci:", subset.locus)
    subset.locus <- paste(subset.locus, collapse = "|")
    tabs <- tabs[grep(subset.locus, tabs)]
  }
  
  ## User-defined exclusion of loci (optional)
  ## ------------------------------------------
  if (!missing(exclude.locus)){
    cat("\nExcluding locus:", exclude.locus)
    exclude.locus <- paste(exclude.locus, collapse = "|")
    tabs <- tabs[-grep(exclude.locus, tabs)]
  }
  
  #############################################
  ##  PART 2: Import alignments
  #############################################
  
  ## Select and read individual alignments
  ## -------------------------------------
  cat("\nReading", length(tabs), "alignments ... ")
  if (missing(subset.species)){
    x <- lapply(tabs, dbReadMSA, x = megProj, blocks = blocks, reliability = reliability)
  } else {
    x <- lapply(tabs, dbReadMSA, x = megProj, taxon = subset.species, regex = FALSE, 
                blocks = blocks, reliability = reliability)
  }
  cat("OK")
  if (any(sapply(x, is.null))) x <- x[!sapply(x, is.null)]
  
  ## Handle blocks
  ## -------------
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
  
  ## Trim tapering ends of alignment
  ## -------------------------------
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
  
  ## Exclude species from alignments that have
  ## less than 'coverage.locus' percent sites
  ## ----------------------------------------
  cat("\nExcluding snippets (<", coverage.locus, "% sites) ... ")
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
  nb.lost <- length(lost)
  if (nb.lost){
    if (nb.lost > 12) lost <- c(head(lost), paste0("[", nb.lost - 12, " sequences]"), tail(lost))
    cat("WARNING:", nb.lost, "species lost:",
        paste("\n-", lost))
  } else {
    cat("OK")
  }
  
  ## Core set of species (optional)
  ## ------------------------------
  if (!missing(core.locus)){
    cat("\nMaking core dataset ... ")
    core.species <- lapply(x[core.locus], function(x) rownames(x))
    core.species <- unique(unlist(core.species))
    subset.alignment <- function(a, s) deleteEmptyCells(a[rownames(a) %in% s, ], quiet = TRUE)
    x <- lapply(x, subset.alignment, s = core.species)
    cat("OK")
  }
  
  ## create partitions
  ## -----------------
  if (!missing(partition)){
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
  if (nrow(p) > 1){
    for (i in 2:nrow(p)){
      p[i, 1] <- p[i - 1, 2] + 1
      p[i, 2] <- p[i, 1] + p[i, 2] -1
    }
  }
  ## partitions in RAxML format:
  p <- paste0("DNA, ", rownames(p), " = ", p[, 1], "-", p[, 2])
  
  ## create SUPERMATRIX
  ## ------------------
  cat("\nConcatenating alignments ... ")
  x <- do.call(cbind.DNAbin, c(x, fill.with.gaps = TRUE))
  cat("OK")
  
  ## exclude species by user decision
  ## --------------------------------
  if (!missing(exclude.species)){
    exclude.species <- intersect(exclude.species, rownames(x))
    nes <- length(exclude.species)
    if (nes > 0){
      cat("\nExcluding ", nes ," (", round(nes/nrow(x), 2), 
          "%) species by user decision ... ", sep = "")
      x <- x[!rownames(x) %in% exclude.species, ]
      cat("OK")
    }
  }
  
  ## Create denser ingroup
  ## ---------------------
  if (best.sampled.congeneric){
    cat("\nKeeping only one best-sampled species per genus ... ")
    percentInformative <- function(z){
      length(which(!z %in% as.raw(c(n = 240, "?" = 2, "-" = 4))))/length(z)
    }
    bsc <- apply(x, 1, percentInformative)
    bsc <- data.frame(species = names(bsc), 
                      genus = strip.spec(names(bsc)),
                      fraction = bsc,
                      stringsAsFactors = FALSE)
    bsc <- bsc[order(bsc$species), ]
    bsc <- split(bsc, f = bsc$genus)
    bsc <- sapply(bsc, function(z) z$species[which.max(z$fraction)])
    x <- x[bsc, ]
    cat("OK")
  }
  
  ## Outgroup
  ## --------
  tax <- dbReadTaxonomy(megProj)
  outgroup <- lapply(megProj@taxon@outgroup, taxdumpChildren,
                     tax = tax, tip.rank = "species")
  outgroup <- do.call(rbind, outgroup)
  outgroup <- intersect(gsub(" ", "_", outgroup$taxon), rownames(x))
  cat("\nNumber of available outgroup species:", length(outgroup))
  
  if (!missing(squeeze.outgroup)){
    cat("\nCreating denser outgroup ... ")
    o <- deleteEmptyCells(x[outgroup, ], quiet = TRUE)
    nn <- as.raw(c(240, 2, 4))
    names(nn) <- c("n", "?", "-")
    nn <- apply(o, 1, function(x, n) length(which((x %in% nn))), n = nn)
    nn <- sort(nn, decreasing = TRUE)
    outgroup <- names(tail(nn, squeeze.outgroup))
    nn <- head(nn, -squeeze.outgroup)
    x <- x[!rownames(x) %in% names(nn), ]
    ## This is a bug, because it changes partitions!
    # x <- deleteEmptyCells(x, quiet = TRUE)
    cat("OK")
  }
  
  ## Make filenames (from here on 'x' does not change any more)
  ## ----------------------------------------------------------
  fn <- paste("data/supermatrix", nrow(x), ngene, ncol(x), sep = "-")
  ext <- c("tre", "phy", "nex", "partitions", "outgroup")
  fns <- paste(fn, ext, sep = ".")
  names(fns) <- ext
  
  ## Assess amount of missing data
  ## -----------------------------
  md <- length(which(x %in% as.raw(c(n = 240, "?" = 2, "-" = 4))))/length(x)
  cat("\nSupermatrix contains ", round(md * 100, 2), 
      "% missing data (including gaps)", sep = "")
  
  ## prepare guide tree 
  ## ------------------
  cat("\nPreparing comprehensive guidetree ... ")
  gt <- comprehensiveGuidetree(megProj, 
                               tip.rank = tip.rank,
                               subset = x)
  cat("OK")
  
  ## write outgroup + partitions files
  ## ---------------------------------
  og <- paste(outgroup, collapse = ",")
  clip <- pipe("pbcopy", "w")
  write(og, file = clip)
  close(clip)
  write(og, fns["outgroup"])
  write(p, fns["partitions"])
  
  ## zip
  ## ---
  cat("\nZipping supermatrix ..\n")
  zip(zipfile = fn, files = fns[c("phy", "partitions", "outgroup")])
  obj <- list(zipname = fn,
              supermatrix = x,
              outgroup = og,
              partitions = p,
              guide.tree = gt)
  
  ## write data as PHY and NEX
  ## -------------------------
  cat("\nWriting supermatrix to files ... ")
  write.tree(gt, fns["tre"])
  write.phy(x, fns["phy"])
  rownames(x) <- gsub("-", "_", rownames(x))
  write.nex(x, fns["nex"])
  cat("done")
  
  obj
}

