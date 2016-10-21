## This code is part of the megaptera package
## © C. Heibl 2013 (last update 2016-09-23)

prune.phylo.rank <- function(phy, tax, rank = "gen"){
  
  ## Expand taxonomy table if column 'species' is missing.
  ## This happens when GenBank was searched for genera.
  ## --------------------------------------------------
  if ( is.null(tax$spec) ){
    spec <- phy$tip.label
    gen <- levels(tax$genus)[tax$genus]
    new.tax <- data.frame()
    for ( i in seq_along(gen) ){
      id <- grep(paste0(gen[i], "_"), spec)
      if ( length(id) > 0 ){
        new.tax <- rbind(new.tax, 
                         data.frame(spec = spec[id], tax[i, ], 
                                    row.names = NULL))
      }
    } # end of FOR-loop
    tax <- new.tax
  }
  
  ## taxonomy
  ## --------
  tax <- tax[which(tax$spec %in% phy$tip.label), ]
  
  rank <- split(tax$spec, tax[rank])
  ## skip incertae sedis '-'
  ic <- grep("^-$", names(rank))
  if ( length(ic) > 0 ) rank <- rank[-ic]
  ntip <- sapply(rank, length)
  
  for ( i in seq_along(rank) ){
    cn <- names(rank)[i]
    id <- which(phy$tip.label %in% rank[[i]])
    cat("\n", i, ": ", cn, " (", length(id), "):", sep = "")
    
    if ( length(id) == 1 ){
      
      ## Case 1: Genus is monotypic
      ## --------------------------
      cat(" monotypic")
      phy$tip.label[id] <- paste(phy$tip.label[id], " - monotypic")
    } else {
      if ( is.monophyletic(phy, id) ) {
        
        ## Case 2: Genus is monophyletic
        ## -----------------------------
        cat(" monophyletic")
        phy$tip.label[id] <- paste(cn, "-", length(id), "sp.")
        phy <- drop.tip(phy, id[-1])
        # TO DO: set edge length to crown group
      }
      else {
        ext <- nonmonophyletic(phy, id)
        ext <- phy$tip.label[ext]
        if ( length(ext)/length(id) < 1 & 
             length(grep("p[.]p[.]|sp[.]", ext)) == 0 ){
          
          ## Case 3: Genus is paraphyletic
          ## (defined as above in IF-clause!)
          ## --------------------------------
          ext <- paste(ext, collapse = ", ")
          cat(" paraphyletic: includes ", ext, sep = "")
          cn <- paste(cn, "incl.", ext)
          phy$tip.label[id] <- paste(cn, "-", length(id), "sp.)")
          phy <- drop.tip(phy, id[-1])
          # TO DO: set edge length to crown group
          
        } else {
          
          ## Case 3: Genus is polyphyletic
          ## -----------------------------
          pp <- proParte(phy, id)
          drop.nodes <- vector()
          for ( j in seq_along(pp) ){
            ppp  <- pp[[j]]
            nt <- length(ppp) ## number of tips in cluster
            if ( nt == 1 ){
              tn <- gsub("^.+_", "", phy$tip.label[ppp])
              tn <- paste(tn, collapse = ", ")
              tn <- paste0("- ", tn, "/", length(id), " sp.")
              phy$tip.label[ppp] <- paste(cn, "p.p.", tn)
            } else {
              tn <- paste0("- ", nt, "/", length(id), " sp.")
              phy$tip.label[ppp] <- paste(cn, "p.p.", tn)
              drop.nodes <- c(drop.nodes, ppp[-1])
            }
          }
          ## tips must be dropped outside of FOR-loop!
          if ( length(drop.nodes) > 0 ) {
            phy <- drop.tip(phy, drop.nodes)
          }
        }
      } 
    }
  }
  return(phy)
}