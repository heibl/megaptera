## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-11-23)

checkSpecLocus <- function(megapteraProj, stage = "sel",
                              plot = TRUE){
  
  ## join taxonomy and locus tables
  ## ------------------------------
  x <- dbReadLocus(megapteraProj)
  spec <- rownames(x)
  
  ## select columns according to 'colname'
  ## -------------------------------------
  stage <- match.arg(stage, c("gb", "sel"))
  stage <- paste(stage, "_", sep = "")
  cols <- grep(stage, names(x))
  if ( length(cols) == 0 ) stop("no locus available")
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
  for ( i in names(pslist) ){
    pslist[[i]] <- sort(rownames(x)[x[, i] == 1])
  }
  
  ## produce barplot
  ## ---------------
  if ( plot ) {
    xx <- t(obj)
    xx[1, ] <- xx[1, ] - xx[2, ]
    opar <- par(no.readonly = TRUE)
    par(mar = c(4, 8, 4, 2))
    taxon <- match.arg(megapteraProj@taxon@tip.rank, 
                       c("genera", "species"))
    df.bar <- barplot(xx, horiz = TRUE, las = 1, xlim = c(0, max(total)),
                      xlab = paste("Number of", taxon),
                      col = c("steelblue1", "orange"), border = NA,
                      main = paste("Number of", taxon, "per locus") #,
                      #                       legend.text = c("species represented by at least two markers", 
                      #                                       "private (only this marker)")
    )
    id <- which(xx[2, ] > 0 )
    if ( length(id) > 0) {
      x <- (2 * xx[1, id] + xx[2, id]) / 2
      y <- df.bar[id]
      v <- xx[2, id]
      text(x, y, v)
    }
    ## plot total number of species as line
    ## makes sense only for more than one locus
    ## ----------------------------------------
    if (  ncol(xx) > 1){
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