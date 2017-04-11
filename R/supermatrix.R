## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-03-08)

# TO DO: marker-wise deletion of species (line 38)

#' @title Retrieval and Concatentation of Multiple Sequence Alignments
#' @description Concatenates multiple sequence alignments of individual loci into a supermatrix.
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @param min.n.seq Numeric, the minimum number of sequences in any alignment that is required to be included in the supermatrix.
#' @param masked Logical, if \code{TRUE} masked (if available) alignments will be concatenated; see \code{\link{stepH}}.
#' @param partition A named list, each element of which contains the names of the loci to group together in one partition. The default is to assign each locus to its own partition. Partitioning of first/second and third nucleotides positions is not yet possible, but might be implemented in the future.
#' @param subset.locus \emph{Currently unused.}
#' @param coverage.locus Numeric between 0 and 1 giving the required minimum coverage of any species in any alignment.
#' @param exclude.locus A vector of mode \code{"character"} giving the names of the loci in the database that will be excluded from the concatenation.
#' @param exclude.species \emph{Currently unused.}
#' @param core.locus A vector of mode \code{"character"} giving the names of the 'core' loci: The resulting supermatrix will only contain species that are contained in the 'core' loci. This option is intended to create denser supermatrices.
#' @param core.species \emph{Currently unused.}
#' @param squeeze.outgroup Numeric, can be given to reduce the number of outgroup species: The function will select the \code{squeeze.outgroup} outgroup species with the best coverage. The idea is to create a more densely sampled outgroup.
#' @return a list with four elements:
#' \item{zipname}{a character string giving the name of the files produced.}
#' \item{supermatrix}{a matrix of class \code{\link{DNAbin}}.}
#' \item{outgroup}{a vector of mode \code{"character"} giving the names of the species used as an outgroup.}
#' \item{partitions}{a vector of mode \code{"character"} giving the partitions of the supermatrix in the same spelling as accepted by RAxML.}
#'
#' In addition, three files are written to the working directory: (1) a NEXUS-formatted file and (2) a PHYLIP-formatted file containing the sequence matrix and (3) a zipped directory containing the PHYLIP-formatted sequence matrix plus output and partitions as separate ASCII-formatted files.
#' @export

supermatrix <- function(megProj, min.n.seq = 3, 
                        masked = TRUE, 
                        partition,
                        subset.locus, coverage.locus = 0.5,
                        exclude.locus, exclude.species,
                        core.locus, core.species,
                        squeeze.outgroup){
  
  ## CHECKS
  ## ------
  if (!inherits(megProj, "megapteraProj")) stop("'megProj' is not of class 'megapteraProj'")
  
  tip.rank <- megProj@taxon@tip.rank
  if (length(min.n.seq) == 1) min.n.seq <- rep(min.n.seq, 2)
  
  ## determine tables that contain more than min.n.seq species
  ## ---------------------------------------------------------
  tabs <- checkBlocks(megProj, plot = FALSE)
  id <- sapply(tabs, function(z) any(z >= min.n.seq[1]))
  tabs <- names(tabs)[id]
  
  ## subset loci to locus
  ## ---------------------
  if (!missing(subset.locus)){
    cat("\n.. subsetting to loci:", subset.locus, " ..")
    subset.locus <- paste(subset.locus, collapse = "|")
    tabs <- tabs[grep(subset.locus, tabs)]
  }
  
  ## exclude loci from concatenation
  ## -------------------------------
  if (!missing(exclude.locus)){
    cat("\n.. excluding locus:", exclude.locus, " ..")
    exclude.locus <- paste(exclude.locus, collapse = "|")
    tabs <- tabs[-grep(exclude.locus, tabs)]
  }
  
  ## select and read individual alignments
  ## -------------------------------------
  cat("\n.. reading", length(tabs), "alignments ..")
  x <- lapply(tabs, dbReadDNA, x = megProj, blocks = "split", masked = masked)

  ## handle blocks
  ## -------------
  block.id <- which(sapply(x, is.list))
  cat("\n..", length(block.id), "alignments are split into blocks ..")
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
    cat("\n.. keeping only blocks >", min.n.seq[2], "sequences ..")
    x[block.id] <- lapply(x[block.id], concatenateBlocks, n = min.n.seq[2])
  }
  
  
  # x <- lapply(x, trimEnds, min.n.seq = min.n.seq)
  
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
  if (length(lost)){
    cat("\n   WARNING:", length(lost), "species lost:",
        paste("\n  -", lost))
  }
  
  ## core set of species
  ## -------------------
  if (!missing(core.locus)){
    cat("\n.. making core dataset ..")
    core.species <- lapply(x[core.locus], function(x) rownames(x))
    core.species <- unique(unlist(core.species))
    subset.alignment <- function(a, s) deleteEmptyCells(a[rownames(a) %in% s, ], quiet = TRUE)
    x <- lapply(x, subset.alignment, s = core.species)
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
  for (i in 2:nrow(p)){
    p[i, 1] <- p[i - 1, 2] + 1
    p[i, 2] <- p[i, 1] + p[i, 2] -1
  }
  ## partitions in RAxML format:
  p <- paste0("DNA, ", rownames(p), " = ", p[, 1], "-", p[, 2])
  
  ## create SUPERMATRIX
  ## ------------------
  cat("\n.. concatenating alignments ..")
  x <- do.call(cbind.DNAbin, c(x, fill.with.gaps = TRUE))
  
  ## exclude species by user decision
  ## --------------------------------
  if (!missing(exclude.species)){
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
  tax <- dbReadTaxonomy(megProj)
  outgroup <- lapply(megProj@taxon@outgroup, taxdumpDaughters,
                     x = tax, tip.rank = "species")
  outgroup <- do.call(rbind, outgroup)
  outgroup <- intersect(gsub(" ", "_", outgroup$taxon), rownames(x))
  cat("\n.. number of available outgroup species:", length(outgroup),  "..")
  
  if (!missing(squeeze.outgroup)){
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
  gt <- comprehensiveGuidetree(megProj, 
                               tip.rank = tip.rank,
                               subset = x)
  
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
  cat("\n.. zipping supermatrix ..\n")
  zip(zipfile = fn, files = fns[c("phy", "partitions", "outgroup")])
  obj <- list(zipname = fn,
              supermatrix = x,
              outgroup = og,
              partitions = p,
              guide.tree = gt)
  
  ## write data as PHY and NEX
  ## -------------------------
  cat("\n.. writing supermatrix to files ..")
  write.tree(gt, fns["tre"])
  write.phy(x, fns["phy"])
  rownames(x) <- gsub("-", "_", rownames(x))
  write.nex(x, fns["nex"])
  
  
  obj
}

