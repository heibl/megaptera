## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2021-03-19)

# TO DO: marker-wise deletion of species (line 38)

#' @title Retrieval and Concatentation of Multiple Sequence Alignments
#' @description Concatenates multiple sequence alignments of individual loci
#'   into a supermatrix.
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @param min.n.seq Numeric, the minimum number of sequences in any alignment
#'   that is required to be included in the supermatrix.
#' @inheritParams pg2DNAbin
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
#' @param locus.coverage Numeric between 0 and 1 giving the required minimum
#'   coverage of any species in any alignment.
#' @param global.coverage Numeric between 0 and 1 giving the required minimum
#'   coverage of any species the concatenated alignment.
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
#' @param protect.outgroup Logical, if \code{TRUE}, the effects of argument
#'   \code{core.locus} and \code{global.coverage} on outgroup taxa will be
#'   ignored.
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
#' @importFrom ape cbind.DNAbin write.tree
#' @importFrom crayon %+% bold cyan magenta red silver
#' @importFrom utils zip
#' @export

supermatrix <- function(megProj, min.n.seq = 3, blocks = "split", 
                        row.confid = 0, col.confid = 0,
                        partition, trim.ends = 0,
                        locus.coverage = 0.0,
                        global.coverage = 0.0,
                        subset.locus = NULL, subset.species = NULL,
                        exclude.locus = NULL, exclude.species = NULL,
                        core.locus = NULL, core.species = NULL,
                        best.sampled.congeneric = FALSE,
                        protect.outgroup = FALSE, squeeze.outgroup = NULL){
  
  ## INITIAL CHECKS + ADJUSTMENTS
  ## ----------------------------
  if (!inherits(megProj, "megapteraProj")) stop("'megProj' is not of class 'megapteraProj'")
  tip.rank <- megProj@taxon@tip.rank
  if (length(min.n.seq) == 1) min.n.seq <- rep(min.n.seq, 2)
  blocks <- match.arg(blocks, c("concatenate", "ignore", "split"))
  if (blocks == "concatenate") blocks <- "split"
  
  ## Get outgroup
  ## ------------
  og <- outgroup(megProj, sep = "_")
  
  x <- selectMSA(megProj = megProj,
                 min.n.seq = min.n.seq,
                 row.confid = row.confid, col.confid = col.confid,
                 trim.ends = trim.ends,
                 locus.coverage = locus.coverage, global.coverage = global.coverage,
                 subset.locus = subset.locus, subset.species = subset.species,
                 exclude.locus = exclude.locus, exclude.species = exclude.species,
                 core.locus = core.locus, core.species = core.species,
                 best.sampled.congeneric = best.sampled.congeneric,
                 protect.outgroup = protect.outgroup, squeeze.outgroup = squeeze.outgroup)
  
  ## Prepare column weights. Note that these must be integers for RAxML
  ## ------------------------------------------------------------------
  cat(silver("Preparing column weights ... "))
  w <- unlist(lapply(x, "attr", which = "cs"))
  w <- round(w / min(w))
  cat(green("OK\n"))
  
  ## Create partitions
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
  cat(silver("Concatenating alignments ... "))
  x <- do.call(cbind.DNAbin, c(x, fill.with.gaps = TRUE))
  cat(green("OK\n"))
  
  ## Exclude species by user decision
  ## --------------------------------
  if (!missing(exclude.species)){
    exclude.species <- intersect(exclude.species, rownames(x))
    nes <- length(exclude.species)
    if (nes > 0){
      cat("Excluding ", nes ," (", round(nes/nrow(x), 2), 
          "%) species by user decision ... ", sep = "")
      x <- x[!rownames(x) %in% exclude.species, ]
      cat(green("OK\n"))
    }
  }
  
  ## Create denser ingroup
  ## ---------------------
  if (best.sampled.congeneric){
    cat("Keeping only one best-sampled species per genus ... ")
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
    cat(green("OK\n"))
  }
  
  ## Detect actual outgroup
  ## ----------------------
  og <- intersect(og, rownames(x))
  cat(silver("Number of available outgroup species: " 
             %+% red$bold(length(og)) %+% "\n"))
  
  ## Make filenames (from here on 'x' does not change any more)
  ## ----------------------------------------------------------
  fn <- paste("data/supermatrix", nrow(x), ngene, ncol(x), sep = "-")
  ext <- c("tre", "phy", "nex", "partitions", "weights", "outgroup")
  fns <- paste(fn, ext, sep = ".")
  names(fns) <- ext
  
  ## Assess amount of missing data
  ## -----------------------------
  md <- length(which(x %in% as.raw(c(n = 240, "?" = 2, "-" = 4))))/length(x)
  cat(silver("Supermatrix (" 
             %+% magenta(bold(nrow(x)) %+% " taxa/" %+% bold(ncol(x)) %+% " characters") 
             %+% ") contains " %+% magenta$bold(paste0(round(md * 100, 2), "%"))
             %+% " missing data (including gaps)\n"))
  
  ## prepare guide tree 
  ## ------------------
  cat("Preparing comprehensive guidetree ... ")
  gt <- comprehensiveGuidetree(megProj, 
                               tip.rank = tip.rank,
                               subset = x)
  cat(green("OK\n"))
  
  ## write outgroup + partitions files
  ## ---------------------------------
  og <- paste(og, collapse = ",")
  clip <- pipe("pbcopy", "w")
  write(og, file = clip)
  close(clip)
  write(og, fns["outgroup"])
  write(p, fns["partitions"])
  
  ## zip
  ## ---
  cat("Zipping supermatrix ...")
  dump <- zip(zipfile = fn, files = fns[c("phy", "partitions", "outgroup")])
  obj <- list(zipname = fn,
              supermatrix = x,
              outgroup = og,
              partitions = p,
              weights = w,
              guide.tree = gt)
  cat(green("OK\n"))
  
  ## write data as PHY and NEX
  ## -------------------------
  cat("\nWriting supermatrix to files ... ")
  write.tree(gt, fns["tre"])
  write.phy(x, fns["phy"])
  rownames(x) <- gsub("-", "_", rownames(x))
  write.nex(x, fns["nex"])
  write(paste(w, collapse = " "), file = fns["weights"]) ## weight file
  cat("done")
  
  obj
}

