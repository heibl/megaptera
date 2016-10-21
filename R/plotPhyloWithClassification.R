## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2016-02-19)

plotPhyloWithClassification <- function(x, phy, og, 
                                        surrogate.tips,
                                        min.support, 
                                        collapse.unsupported = TRUE,
                                        ignore.tag,
                                        edges, 
                                        x.lim,
                                        pdf.height = "auto",
                                        pdf.width = 20,
                                        title, file){
  
  ## ignore tag
  if ( !missing(ignore.tag) ){
    id <- grep(ignore.tag, phy$tip.label, fixed = TRUE)
    phy$tip.label[id] <- tagged  <- gsub(paste("_", ignore.tag, sep = ""), "", 
                                         phy$tip.label[id],
                                         fixed = TRUE)
  }
  
  ## modify tree
  ## -----------
  if ( !missing(og) ) {
    if ( is.monophyletic(phy, og) ) {
      phy <- root(phy, og)
      phy <- drop.tip(phy, og)
    } else {
      phy <- drop.tip(phy, og)
    }
  }
  phy <- ladderize(phy)
  phy <- fixNodes(phy)
  
  ## append asterisk to species that are
  ## represented by a surrogate
  ## --------------------------
  tax <- dbReadTaxonomy(x, subset = phy)
  s <- tax$spec[grep("surrogate", tax$tag)]
  if ( length(s) > 0 ){
    id.surrogate <- which(phy$tip.label %in% s)
  } else {
    id.surrogate <- NULL
  }
  id.mtra <- grep("[+]", phy$tip.label)
  
  phy$tip.label <- gsub("_[+|*]", "", phy$tip.label)
  tax <- dbReadTaxonomy(x, subset = phy) ## needs to be updated!
  
  cat(phy$tip.label)
  
  ## monophyly of genera
  ## -------------------
  if ( x@taxon@tip.rank == "spec" ){
    gen <- split(tax$spec, tax$gen)
    id <- sapply(gen, is.monophyletic, phy = phy)
    if ( missing(edges) ){
      ecol <- edge.color(phy, gen[id], col = "#2166ac")
    } else {
      ecol <- edge.color(phy, edges, col = "orange")
    }
    
    tcol <- tip.color(phy, gen[!id], regex = TRUE, 
                      col = rep("#b2182b", length(which(!id))))
    mp.gen <- noi(phy, gen[id])
    n.spec <- sapply(gen[id], length)
    gen.names <- names(mp.gen)
    gen.names[n.spec == 1] <- ""
    
    ## for monophyletic genera with more than 1 species,
    ## display only epithet:
    g <- gen.names[gen.names != ""]
    g <- phy$tip.label %in% unlist(gen[g])
    
  } else {
    ecol <- tcol <- "black"
  }
  
  ## colors for families
  ## -----------------
  fam <- split(tax[, x@taxon@tip.rank], tax$fam)
  id <- sapply(fam, is.monophyletic, phy = phy)
  mp.fam <- noi(phy, fam[id])
  
  if ( x@taxon@tip.rank == "spec" ){
    phy$tip.label[g] <- gsub("^.+_", "", phy$tip.label[g]) ## cf. line 44
  }
  
  ## representation by markers
  ## -------------------------
  #   mm <- check.Coverage(conn, pool.markers = FALSE)
  #   mm <- mm[match(phy$tip.label, rownames(mm)), -(1:2)]
  #   names(mm) <- gsub("X_|_sel", "", names(mm))
  #   # mm <- mm[, -which(colSums(mm) == 0)]
  #   mm <- mm + 1
  #   bubble <- function(x, phy){
  #     for ( i in seq_along(names(x)) )
  #       append2tips(phy, align = TRUE, grid = TRUE, 
  #                   pch = 19, col = c("grey50", "red")[x[, i]], legend = names(x)[i])
  #   }
  #   
  
  ## reappend tags
  if ( length(id.surrogate) > 0 )
    phy$tip.label[id.surrogate] <- paste(phy$tip.label[id.surrogate], 
                                         "*", sep = "_")
  if ( length(id.mtra) > 0 )
    phy$tip.label[id.mtra] <- paste(phy$tip.label[id.mtra], 
                                    "+", sep = "_")
  if ( !missing(ignore.tag) ){
    id <- phy$tip.label %in% tagged
    phy$tip.label[id] <- paste(phy$tip.label[id], ignore.tag, sep = "_")
  }
  
  
  if ( collapse.unsupported )
    phy <- collapseUnsupportedEdges(phy, cutoff = min.support)
  
  ## open PDF device
  ## ---------------
  if ( pdf.height == "auto" ) pdf.height <- Ntip(phy)/5
  pdf(file, 
      # paper = "a4", 
      height = pdf.height, width = pdf.width)
  plot(phy, label.offset = .001, cex = .75, 
       tip.color = "white", edge.color = NA,
       x.lim = x.lim
  )
  box.clades(phy, mp.fam, col = "#d1e5f0", text = names(mp.fam),
             align = "all")
  if ( x@taxon@tip.rank == "spec" ){
    box.clades(phy, mp.gen, col = "#67a9cf", text = gen.names)
  }
  plot.phylo.upon(phy, label.offset = .001, cex = .75, 
                  tip.color = tcol, edge.color = ecol,
                  x.lim = x.lim
  )
  #bubble(mm, phy) # representation by markers
  #append2tips(phy, align = phyUE, grid = TRUE, pch = 19)
  
  if ( !missing(min.support) ) 
    node.support(phy$node.label, cutoff = min.support, cex = .6)
  add.scale.bar(length = .05, lwd = 1.5)
  if ( !missing(title) ) title(main = title)
  # legend("topleft", inset = c(.2, .35), legend = names(tax), 
  #        fill = clrs, ncol = 3)
  dev.off()
  invisible(lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv))
}



