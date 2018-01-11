## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-11-29)

#' @title Barplot of Species Numbers per Locus
#' @description Produce a barplot showing the number of species found for each
#'   locus. Thereby the number of species either refers to the number of species
#'   found on GenBank, the number of species selected according to identity and
#'   coverage, or the number of species included in the final alignment.
#'   Additionally, black lines to right of the bars indicate the cumulative
#'   species numbers for each locus if this locus was to be concatenated with
#'   all other loci having more species.
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @param stage A vector of mode \code{character}; valid are \code{"gb"},
#'   \code{"sel"}, or \code{"blocks"} referring to the number of species (1)
#'   found on GenBank, (2) selected according to coverage and identity (see
#'   \code{\link{locus}}); or included in the final alignment, respectively.
#' @param subset A vector of mode \code{"character"}, that can be used to choose
#'   a subset of the total taxa available.
#' @param plot Logical, indicating if the barplot should be produced on the
#'   current graphical device or not.
#' @return A list, but mainly called for its side effect, the plotting of a
#'   barplot.
#' @seealso \code{\link{checkMissingSpec}} to create a table of species that
#'   were missed/lost during consecutive steps of the pipeline.
#' @importFrom graphics barplot lines par text
#' @export

checkSpecLocus <- function(megProj, stage = "sel", subset,
                           plot = TRUE){
  
  ## join taxonomy and locus tables
  ## ------------------------------
  x <- dbReadLocus(megProj, subset = subset)
  spec <- rownames(x)
  
  ## select columns according to 'colname'
  ## -------------------------------------
  stage <- match.arg(stage, c("gb", "sel"))
  stage <- paste(stage, "_", sep = "")
  cols <- grep(stage, names(x))
  if (!length(cols)) stop("no locus available")
  x <- x[, cols, drop = FALSE]
  
  ## some cosmetics on locus names
  ## -----------------------------
  colnames(x) <- gsub(stage, "", colnames(x))
  colnames(x) <- gsub("^_", "", colnames(x))
  colnames(x) <- gsub("([[:digit:]])(_)([[:digit:]])", "\\1,\\3", colnames(x))
  colnames(x) <- gsub("([[:alpha:]])(_)([[:alpha:]])", "\\1 \\3", colnames(x))
  
  ## convert to binary matrix
  ## ------------------------
  bincov <- function(x){
    x[is.na(x)] <- 0
    x[grep("excluded", x)] <- 0
    x[grep("selected", x)] <- 1
    x[x > 1] <- 1
    as.numeric(x)
  }
  x <- apply(x, 2, bincov)
  rownames(x) <- spec
  
  ## delete species with any locus
  ## -----------------------------
  x <- x[rowSums(x) != 0, , drop = FALSE]
  
  ## number of species per locus
  ## ---------------------------
  sfreq <- sort(colSums(x), decreasing = TRUE)
  x <- x[, match(names(sfreq), colnames(x)), drop = FALSE]
  total <- vector()
  for (i in 1:ncol(x)){
    total <- c(total, length(which(rowSums(x[, 1:i, drop = FALSE]) != 0)))
  }
  
  ## number of private species per locus
  ## -----------------------------------
  mfreq <- rowSums(x)
  private <- which(mfreq == 1)
  x <- x[private, , drop = FALSE]
  pfreq <- colSums(x)
  obj <- cbind(Ntotal = sfreq, Nprivate = pfreq[match(names(sfreq), names(pfreq))])
  
  ## number of loci per species
  ## --------------------------
  pslist <- as.list(names(pfreq)[pfreq > 0])
  names(pslist) <- as.vector(pslist)
  for (i in names(pslist)){
    pslist[[i]] <- sort(rownames(x)[x[, i] == 1])
  }
  
  ## produce barplot
  ## ---------------
  if (plot) {
    xx <- t(obj)
    xx[1, ] <- xx[1, ] - xx[2, ]
    opar <- par(no.readonly = TRUE)
    par(mar = c(4, 8, 4, 2))
    taxon <- ifelse(megProj@taxon@tip.rank == "genus", "genera", "species") 
    df.bar <- barplot(xx, horiz = TRUE, las = 1, xlim = c(0, max(total)),
                      xlab = paste("Number of", taxon),
                      col = c("steelblue1", "orange"), border = NA,
                      main = paste("Number of", taxon, "per locus") #,
                      #                       legend.text = c("species represented by at least two markers", 
                      #                                       "private (only this marker)")
    )
    id <- which(xx[2, ] > 0)
    if (length(id)) {
      x <- (2 * xx[1, id] + xx[2, id]) / 2
      y <- df.bar[id]
      v <- xx[2, id]
      text(x, y, v)
    }
    ## plot total number of species as line
    ## makes sense only for more than one locus
    ## ----------------------------------------
    if (ncol(xx) > 1){
      xyt <- data.frame(total, total, df.bar, df.bar)
      xyt[, 3] <- xyt[, 3] - 0.5 * (df.bar[2] - df.bar[1])
      xyt[, 4] <- xyt[, 4] + 0.5 * (df.bar[2] - df.bar[1])
      for ( i in 2:nrow(xyt) ) lines(xyt[i, 1:2], xyt[i, 3:4], 
                                     lwd = 1, col = "black", lty = 1, xpd = NA)
      text(total[-1], df.bar[-1], total[-1], pos = 2)
      text(total[1], df.bar[1], total[1], pos = 4)
      par(opar)
    }
  }
  
  ## return results
  ## --------------
  list(specPerMarker = obj, privateSpec = pslist)
}