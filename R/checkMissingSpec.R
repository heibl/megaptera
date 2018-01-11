## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-01-11)

#' @title Check for Missing Species
#' @description Create a table of species that failed to pass steps in the
#'   pipeline.
#' @param x An object of class \code{\linkS4class{dbPars}} or
#'   \code{\linkS4class{megapteraProj}}.
#' @param provenance A character string as used in the column \code{provenance}
#'   the \code{acc_<locus>} tables in the postgreSQL database. Depending on your
#'   actual pipeline setup these may be \code{"ncbi"}, \code{"ncbi annotated"},
#'   or \code{"bold"} as well as any user-defined strings.
#' @param tag A character string as used in the column \code{tag} in the
#'   \code{taxonomy} table in the postgreSQL database. Depending on your choices
#'   for parameter values in \code{\link{taxon}} or \code{\link{taxonGuidetree}}
#'   these may be \code{"ingroup (NCBI)"}, \code{"extended ingroup (NCBI)"},
#'   \code{"outgroup (NCBI)"}, or \code{"extended outgroup (NCBI)"} as well as
#'   any user-defined strings.
#' @return A data frame with two columns: 
#' \describe{ 
#'   \item{spec}{Names of missing species} 
#'   \item{status}{Why species is missing; might be (1)
#'   \code{not.on.genbank} (species for which no sequences were found on
#'   GenBank), (2) \code{not.selected} (species that were not selected due to
#'   below-threshold identity or coverage), or (3) \code{not.aligned} (species
#'   that were lost during the alignment process; happens rarely)}
#'   }
#' @seealso \code{\link{checkSpecLocus}} to create a barplot of the species number per locus.


checkMissingSpec <- function(x, provenance = ".", tag = "ingroup"){  
  
  ## get locus table:
  cat("tag        :", tag, "\nprovenance :", provenance)
  lcs <- dbReadLocus(x, provenance = provenance, tag = tag)

  
  ## species/genera and their numbers:
  spec <- total <- rownames(lcs)
  nspec <- length(spec)
  gen <- unique(strip.spec(spec))
  ngen <- length(gen)
  
  ## table entries must be coerced to {0, 1}:
  make.binary <- function(y, sel0, sel1){
    y[is.na(y)] <- 0 
    if ( !missing(sel0) ) y[grep(sel0, y)] <- 0
    if ( !missing(sel1) ) y[grep(sel1, y)] <- 1
    y <- as.numeric(y)
    y[y > 1] <- 1
    y
  }
#   lcs <- apply(lcs, 2, make.binary)
#   rownames(lcs) <- spec ## make.binary strips rownames!
  
  ## Not found on GenBank
  ## --------------------
  cat("\n\nno sequences available:")
  ## species:
  id <- c(TRUE, FALSE)
  lcs[, id] <- apply(lcs[, id], 2, make.binary)
  rownames(lcs) <- spec
  gb <- rowSums(lcs[, id, drop = FALSE])
  gb <- sort(names(gb)[gb == 0])
  ## genera:
  present.species <- setdiff(spec, gb)
  present.genera <- strip.spec(present.species)
  missing.genera <- strip.spec((gb))
  gb.gen <- setdiff(missing.genera, present.genera)
  ## screen output:
  gb.freq <- format(c(length(gb), length(gb.gen)))
  gb.frac <- c(length(gb) / nspec, length(gb.gen) / ngen)
  gb.frac <- round(gb.frac, 2)
  cat("\n .. ", gb.freq[1], " of ", nspec, " species (", 
      gb.frac[1], "%)", 
      "\n .. ", gb.freq[2], " of ", ngen, " genera (", 
      gb.frac[2], "%)", sep = "")
  nspec <- length(present.species)
  ngen <- length(unique(present.genera))
  lcs1 <- lcs2 <- lcs[rownames(lcs) %in% present.species, ]
  spec <- rownames(lcs1)
  
  ## found on GenBank, but not selected
  ## in stepF due to identity or coverage
  ## ------------------------------------
  cat("\nnot selected:")
  ## species:
  id <- c(FALSE, TRUE)
  lcs1[, id] <- apply(lcs1[, id], 2, make.binary, 
                     sel1 = "selected|aligned|masked|excluded")
  rownames(lcs1) <- spec
  sel <- rowSums(lcs1[, id, drop = FALSE])
  sel <- sort(names(sel[sel == 0]))
  ## genera:
  present.species <- setdiff(spec, sel)
  present.genera <- strip.spec(present.species)
  missing.genera <- strip.spec((sel))
  sel.gen <- setdiff(missing.genera, present.genera)
  ## screen output:
  sel.freq <- format(c(length(sel), length(sel.gen)))
  sel.frac <- c(length(sel) / nspec, length(sel.gen) / ngen)
  sel.frac <- round(sel.frac, 2)
  cat("\n .. ", sel.freq[1], " of ", nspec, " species (", 
      sel.frac[1], "%)", 
      "\n .. ", sel.freq[2], " of ", ngen, " genera (", 
      sel.frac[2], "%)", sep = "")
  nspec <- nspec - length(sel)
  ngen <- ngen - length(unique(missing.genera))
  
  ## dropped in the alignment process
  ## --------------------------------
  ## species:
  cat("\nexcluded by user:")
  id <- c(FALSE, TRUE)
  lcs2[, id] <- apply(lcs2[, id], 2, make.binary, sel0 = "excluded",
                      sel1 = "selected|aligned|masked")
  rownames(lcs2) <- spec
  excluded <- rowSums(lcs2[, id, drop = FALSE])
  excluded <- names(excluded[excluded == 0])
  ## genera:
  present.species <- setdiff(spec, excluded)
  present.genera <- strip.spec(present.species)
  missing.genera <- strip.spec((excluded))
  excluded.gen <- setdiff(missing.genera, present.genera)
  ## screen output:
  excluded.freq <- format(c(length(excluded), length(excluded.gen)))
  excluded.frac <- c(length(excluded) / nspec, length(excluded.gen) / ngen)
  excluded.frac <- round(excluded.frac, 2)
  cat("\n .. ", excluded.freq[1], " of ", nspec, " species (", 
      excluded.frac[1], "%)", 
      "\n .. ", excluded.freq[2], " of ", ngen, " genera (", 
      excluded.frac[2], "%)", sep = "")
  
  ## concataneta missing species and genera
  ## --------------------------------------
  spec <- c(gb, sel, excluded)
  gen <- unlist(unique(c(gb.gen, sel.gen, excluded.gen)))
  
  ## genera with only two species = "sisters"
  ## and which have only one species in the database
  ## -----------------------------------------------
  sister <- table(strip.spec(total))
  sister <- names(sister)[sister == 2]
  one.missing <- table(strip.spec(spec))
  one.missing <- names(one.missing)[one.missing == 1]
  sister.one.missing <- intersect(sister, one.missing)
  
  ## prepare summmary object for output
  ## ----------------------------------
  obj <- data.frame(spec = spec,
                    status = c(rep("not on GenBank", length(gb)),
                               rep("not selected", length(sel)),
                               rep("excluded by user", length(excluded))),
                    genus = "uncomplete",
                    stringsAsFactors = FALSE)
  id <- grep(paste(gen, collapse = "|"), spec)
  obj$genus[id] <- "entirelyMissing"
  id <- grep(paste(sister.one.missing, collapse = "|"), spec)
  obj$genus[id] <- "oneOfTwo"
  rownames(obj) <- NULL
  obj
}
