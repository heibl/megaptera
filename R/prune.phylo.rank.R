## This code is part of the megaptera package
## © C. Heibl 2013 (last update 2017-11-13)

#' @title Prune Phylogenies to Higher Ranks
#' @description Prune tips of a certain taxonomic rank from an object of class
#'   \code{\link[ape]{phylo}} and obtain a new \code{phylo} object whose tips
#'   are of a higher rank.
#' @param phy An object of class \code{\link[ape]{phylo}}.
#' @param tax A data frame containing taxonomic information for the tip labels
#'   in \code{phy}.
#' @param rank A character string giving the name of a column (= taxonomic rank)
#'   in \code{tax} to which \code{phy} will be pruned.
#' @return An object of class \code{\link[ape]{phylo}}.
#' @seealso \code{\link{addTips}} and \code{\link{addSingleTip}} to add terminal nodes to a phylogeny.
#' @importFrom ape drop.tip is.monophyletic
#' @export

prune.phylo.rank <- function(phy, tax, rank = "genus"){
  
  ## Expand taxonomy table if column 'species' is missing.
  ## This happens when GenBank was searched for genera.
  ## --------------------------------------------------
  if ( is.null(tax$spec) ){
    spec <- phy$tip.label
    gen <- levels(tax$genus)[tax$genus]
    new.tax <- data.frame()
    for (i in seq_along(gen)){
      id <- grep(paste0(gen[i], "_"), spec)
      if (length(id)){
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
  if (length(ic)) rank <- rank[-ic]
  ntip <- sapply(rank, length)
  
  for (i in seq_along(rank)){
    cn <- names(rank)[i]
    id <- which(phy$tip.label %in% rank[[i]])
    cat("\n", i, ": ", cn, " (", length(id), "):", sep = "")
    
    if (length(id) == 1){
      
      ## Case 1: Genus is monotypic
      ## --------------------------
      cat(" monotypic")
      phy$tip.label[id] <- paste(phy$tip.label[id], " - monotypic")
    } else {
      if (is.monophyletic(phy, id)) {
        
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
        if (length(ext)/length(id) < 1 & 
             length(grep("p[.]p[.]|sp[.]", ext)) == 0){
          
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
          
          ## Case 4: Genus is polyphyletic
          ## -----------------------------
          pp <- proParte(phy, id)
          drop.nodes <- vector()
          for (j in seq_along(pp)){
            ppp  <- pp[[j]]
            nt <- length(ppp) ## number of tips in cluster
            if (nt == 1){
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
          if (length(drop.nodes)) {
            phy <- drop.tip(phy, drop.nodes)
          }
        }
      } 
    }
  }
  return(phy)
}