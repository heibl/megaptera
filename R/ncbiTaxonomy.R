## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-12-17)

ncbiTaxonomy <- function(taxon, kingdom, extend = FALSE, trim, 
                         megapteraProj){
  
 kingdom <- match.arg(kingdom, c("Fungi", "Metazoa", "Viridiplantae"))
  
  taxon.untouched <- taxon
  species.list <- unique(is.Linnean(unlist(taxon)))
  if ( length(species.list) > 1 ) 
    stop("names of species and higher taxa", 
         " must not be mixed")
  
  ##
  if ( species.list ){
    accepted <- sapply(taxon, head, 1) # get first element
    syno <- unlist(sapply(taxon, tail, -1))
    accepted2syno <- taxon[sapply(taxon, length) > 1] 
    acceptedsyno <- union(accepted, syno)
    acceptedsyno <- data.frame(
      genus = strip.spec(acceptedsyno), 
      species = acceptedsyno,
      stringsAsFactors = FALSE)
    taxon <- unique(acceptedsyno$genus)
    taxon <- sort(taxon) ## makes debugging easier
  }
  
  ## download in batches of 50
  ## -------------------------
  id <- seq(from = 1, to = length(taxon), by = 50)
  id <- data.frame(from = id, to = c(id[-1] - 1, length(taxon)))
  z <- list()
  for ( i in 1:nrow(id) ){
    j <- seq(from = id$from[i], to = id$to[i])
    zz <- ncbiLineage(taxon = taxon[j],
                      kingdom = kingdom,
                      megapteraProj = megapteraProj)
    if ( length(zz) == 0 ){
      write(taxon[j], file = "ncbiTaxonomy-missing.txt", 
            append = TRUE)
    } else {
      z <- c(z, zz)
    }
  }
  
  ## returned list is empty
  if ( length(z) == 0 ) {
    return(NULL)
  }
  
  ## delete ranks that are incertae sedis
  ## ------------------------------------
  dis <- function(z){
    id <- grep("incertae sedis", z$name)
    if ( length(id) > 0 ){
      z <- z[-id, ]
    }
    z
  }
  z <- lapply(z, dis) 
  
  ## check ranks and add columns if necessary
  ## ----------------------------------------
  ranks <- sortRanks(z)
  z <- lapply(z, addRanks, ranks)
  z <- lapply(z, function(z) z$name)
  z <- do.call(rbind, z)
  colnames(z) <- ranks
  z <- data.frame(z, stringsAsFactors = FALSE)
  
  ## intersect genus-level taxonomy with species list
  ## (necessary, because taxonomic classification was searched
  ## via generic names)
  ## ------------------
  if ( species.list ){
    if ( extend ){
      acceptedsyno <- rbind(acceptedsyno,
                            z[, c("genus", "species")])
      acceptedsyno <- unique(acceptedsyno)
    }
    z <- unique(z[, 1:which(names(z) == "genus")])
    z <- data.frame(z[match(acceptedsyno$genus, z$genus), ], 
                    species = acceptedsyno$species,
                    syno = acceptedsyno$species %in% syno,
                    synonym = "-",
                    stringsAsFactors = FALSE)
    
    ## handle synonyms
    ## ---------------
    if ( length(accepted2syno) > 0 ){
      ## do accepted species with synonyms ahve a taxonomy?
      a <- sapply(accepted2syno, head, 1)
      a <- z[!is.na(z$genus) & z$species %in% a, "species"]
      ## if so, append their synonyms
      if ( length(a) > 0 ){
        as <- function(w) c(head(w, 1), paste(tail(w, -1), collapse = "|"))
        as <- do.call(rbind, lapply(accepted2syno, as))
        z$synonym[match(as[, 1], z$species)] <- as[, 2]
      }
     
      ## do synonyms have a taxonomy?
      s <- unlist(sapply(accepted2syno, tail, -1))
      s <- z[!is.na(z$genus) & z$species %in% s, "species"]
      ## if so: implement
      if ( length(s) > 0){
        stop("implement me!")
      }
    }
    
    ## identify species that are not listed in the NCBI
    ## taxonomy and check if information is available
    ## for their synonym
    ## ----------------------------------------------
    no.info.accepted <- is.na(z$genus) & !z$syno
    if ( any(no.info.accepted) ){
      no.info.accepted <- sort(z$species[no.info.accepted])
      info.syno <- z$species[!is.na(z$genus) & z$syno]
      ## niaws: which accepted species that were not found
      ## have synonyms
      niaws <- sapply(no.info.accepted, grep, x = accepted2syno)
      niaws <- unlist(niaws)
      if ( length(niaws) > 1 ){
        for ( i in niaws ){
          # s: informative(!) synonyms
          s <- intersect(info.syno,
                         tail(accepted2syno[[i]], -1))
          # a: corresonding accepted name without info
          a <- head(accepted2syno[[i]], 1)
          z <- z[!z$species %in% a, ]
          j <- z$species %in% s
          if ( any(j) ){
            z$genus[j] <- strip.spec(a)
            z$species[j] <- a
            z$syno[j] <- FALSE
            z$synonym[j] <- paste(s, collapse = "|")
          } # end of IF-clause
        } # end of FOR-loop
      } # end of IF-clause
    }
    
    ##  delete synonyms and column 'syno'
    z <- z[!z$syno, names(z) != "syno"]
    ## delete species without taxonomy
    z <- z[!is.na(z$genus), ]

    
    ## identify species that are not listed in the NCBI
    ## taxonomy, write them to file and issue warning
    ## ----------------------------------------------
    no.info <- setdiff(sapply(taxon.untouched, head, n = 1),
                       z$species)
    if ( length(no.info) > 0 ){
      write(sort(no.info), 
            file = "ncbiTaxonomy-missing.txt")
      warning("no classification found for ", 
              length(no.info), " species:\n ", 
              paste(sort(no.info), collapse = "\n "),
              "\n\n(this list is written to file 'ncbiTaxonomy-missing.txt')", 
              sep = " ", .call = FALSE)
    }
  }
  cat("\n.. number of species found:", nrow(z), "..")
  
  ## delete internal, noninformative "no rank"-columns
  ## -------------------------------------------------
  id <- apply(z, 2, unique)
  id <- which(id == "-")
  if ( length(id) > 0 ){
    z <- z[, -id]
  }
  
  ## trimming taxonomy to required lower ranks
  ## -----------------------------------------
  if ( !missing(trim) ){
    trim <- match.arg(trim, c("auto", colnames(z)))
    if ( trim == "auto" ){
      id <- max(which(sapply(apply(z, 2, unique), length) == 1))
    }
    if ( trim %in% colnames(z) ){
      id <- which(colnames(z) %in% trim)
    }
    z <- z[, id:ncol(z)]
  } 
  
  ## enforce binomials; sometimes authorities
  ## and years are added to species names
  ## e.g. 'Gaussia princeps H.Wendl. 1865'
  ## -------------------------------------
  z$species <- strip.infraspec(z$species)
  
  ## order taxonomy
  ## --------------
  for ( i in rev(colnames(z)) ){
    z <- z[order(z[, i]), ]
  }
  z
}