## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-12-09)

fixTaxonomy <- function(tax, auto = FALSE, ignore = c("synonym", "tag"),
                        quiet = TRUE){
  
  ## enforce SQL-compatible attribute names
  tax <- sqlTaxonomyHeader(tax)
  
  ## save ignored columns 'ic', 
  ## i.e. columns without taxonomic meaning
  ic <- names(tax) %in% ignore
  if ( any(ic) ){
    tax.ic <- tax[, ic, drop = FALSE]
    rownames(tax.ic) <- tax$spec
    tax <- tax[, !ic]
  }
  
  ## enforce screen dialog in interactive mode
  if ( !auto ) quiet <- FALSE
  
  ## adjust order: rightmost column must be the highest rank
  ## -------------------------------------------------------
  if ( all(c("spec", "gen") %in% names(tax)) ){
    if ( which(names(tax) %in% "spec") > which(names(tax) %in% "gen") )
      tax <- tax[, ncol(tax):1]
  } else {
    if ( !all(c("fam", "gen") %in% names(tax)) ) stop("no columns spec, gen, fam")
    if ( which(names(tax) %in% "gen") > which(names(tax) %in% "fam") )
      tax <- tax[, ncol(tax):1]
  }
  
  ## sometimes no genus name is returned by NCBI
  ## -------------------------------------------
  id <- tax$gen == "-"
  if ( any(id) ){
    gen <- strip.spec(tax$spec[id])
    levels(tax$gen) <- union(levels(tax$gen), gen)
    tax$gen[id] <- gen
  }
  
  for ( i in (ncol(tax) - 1):1 ){
    if ( !quiet ) cat("Taxonomic rank:", names(tax)[i])
    xx <- unique(tax[, ncol(tax):i])
    d <- xx[, ncol(xx)]
    d <- d[!d %in% c("-", "incertae sedis")]
    dd <- duplicated(d)
    if ( any(dd) ){
      ## identify and loop over duplicate taxa
      d <- unique(d[dd])
      for ( j in d){
        if ( !quiet )  cat("\nProblem found in: rank '", names(tax)[i], "', taxon '", j, "'", sep = "") 
        ## identify ambiguous higher rank and loop over them
        y <- xx[xx[, ncol(xx)] == j, ]
        y <- apply(y, 2, unique)
        y <- y[sapply(y, length) > 1]
        if ( !quiet )  cat("\n", length(y), "higher ranks are ambiguous")
        for ( k in names(y) ){
          if ( !quiet )  cat("\nRank '", k ,"': please choose one of ...", sep = "")
          z <- y[[k]]
          if ( !quiet )  cat(paste("\n (", 1:length(z), ") '", z, "'", sep = ""))
          if ( auto ){
            ind <- sample(seq_along(z), 1) ## random
          } else {
            ind <- as.numeric(readline()) ## interactive
          }
          if ( !quiet )  cat(" your choice: '", z[ind], "'\n", sep = "")
          tax[tax[, i] == j, k] <- z[ind]
        } # end of FOR-loop over k
      } # end of FOR-loop over j
    } else {
      if ( !quiet ) cat(" ... OK\n")
      d <- NULL
    }
  } # end of FOR-loop over i
  tax <- unique(tax[, ncol(tax):1])
  
  ## restore ignored columns
  ## -----------------------
  if ( any(ic) ){
    tax <- data.frame(
      tax,
      tax.ic[match(rownames(tax.ic), tax$spec), ])
    names(tax)[(ncol(tax) + 1) - (length(ignore):1)] <- ignore
  }
  tax
}

