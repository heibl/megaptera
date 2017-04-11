## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2017-02-20)

#' @title Plot Phylogentic Tree with Clade Highlighting
#' @description Plots a phylogenetic tree together with color highlighting of
#'   the underlying taxonomic classification.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param phy An object of class \code{\link{phylo}}.
#' @param og A vector of mode \code{"character"}.
#' @param surrogate.tips
#' @param min.support
#' @param collapse.unsupported Logical
#' @param ignore.tag \emph{Currently unused.}
#' @param edges
#' @param higher.rank
#' @param pdf.height
#' @param pdf.with
#' @param title
#' @param file
#' @importFrom grDevices dev.off pdf
#' @export

plotPhyloWithClassification <- function(x, phy, og, 
                                        surrogate.tips,
                                        min.support, 
                                        collapse.unsupported = TRUE,
                                        ignore.tag,
                                        edges, 
                                        higher.rank,
                                        pdf.height = "auto",
                                        pdf.width = 20,
                                        title, file){
  
  ## ignore tag
  ## ----------
  if (!missing(ignore.tag)){
    id <- grep(ignore.tag, phy$tip.label, fixed = TRUE)
    phy$tip.label[id] <- tagged  <- gsub(paste("_", ignore.tag, sep = ""), "", 
                                         phy$tip.label[id],
                                         fixed = TRUE)
  }
  
  ## modify tree
  ## -----------
  if (!missing(og)) {
    if (is.monophyletic(phy, og)) {
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
  tax$taxon <- gsub(" ", "_", tax$taxon)
  # s <- tax$spec[grep("surrogate", tax$tag)]
  # if (length(s)){
  #   id.surrogate <- which(phy$tip.label %in% s)
  # } else {
  #   id.surrogate <- NULL
  # }
  # id.mtra <- grep("[+]", phy$tip.label)
  # 
  # phy$tip.label <- gsub("_[+|*]", "", phy$tip.label)
  # tax <- dbReadTaxonomy(x, subset = phy) ## needs to be updated!
  
  # cat(phy$tip.label)
  
  ## monophyly of genera
  ## -------------------
  if (x@taxon@tip.rank == "spec"){
    gen <- data.frame(gen = strip.spec(phy$tip.label),
                      spec = phy$tip.label,
                      stringsAsFactors = FALSE)
    gen <- split(gen$spec, gen$gen)
    id <- sapply(gen, is.monophyletic, phy = phy)
    if (missing(edges)){
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
  
  ## colors for higher ranks
  ## -----------------------
  if (!missing(higher.rank)){
    hr <- lapply(phy$tip.label, taxdumpLineage, x = tax)
    hr <- sapply(hr, function(z) z$taxon[z$rank == higher.rank])
    hr[!sapply(hr, length)] <- "incertae sedis"
    hr <- data.frame(hr = unlist(hr), spec = phy$tip.label, stringsAsFactors = FALSE)
    hr <- split(hr$spec, hr$hr)
    id <- sapply(hr, is.monophyletic, phy = phy)
    mp.hr <- noi(phy, hr[id])
  }
  
  
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
  # if ( length(id.surrogate) > 0 )
  #   phy$tip.label[id.surrogate] <- paste(phy$tip.label[id.surrogate], 
  #                                        "*", sep = "_")
  # if ( length(id.mtra) > 0 )
  #   phy$tip.label[id.mtra] <- paste(phy$tip.label[id.mtra], 
  #                                   "+", sep = "_")
  # if ( !missing(ignore.tag) ){
  #   id <- phy$tip.label %in% tagged
  #   phy$tip.label[id] <- paste(phy$tip.label[id], ignore.tag, sep = "_")
  # }
  
  
  if (collapse.unsupported)
    phy <- collapseUnsupportedEdges(phy, cutoff = min.support)
  
  ## open PDF device
  ## ---------------
  if (pdf.height == "auto") pdf.height <- Ntip(phy)/5
  xlim <- max(tipHeights(phy)) * 1.05
  pdf(file, 
      # paper = "a4", 
      height = pdf.height, width = pdf.width)
  plot(phy, label.offset = .001, cex = .75, 
       tip.color = "white", edge.color = NA,
       x.lim = xlim
  )
  ## highlight higher rank: color topology
  ## -------------------------------------
  if (!missing(higher.rank)) {
    fillClades(phy, mp.hr, col = "lightsteelblue2")
  }
  if ( x@taxon@tip.rank == "spec" ){
    box.clades(phy, mp.gen, col = "sandybrown", text = gen.names)
  }
  plotPhyloUpon(phy, label.offset = .001, cex = .75, 
                  tip.color = tcol, edge.color = ecol,
                  x.lim = xlim
  )
  ## highlight higher rank: annotate edges
  ## -------------------------------------
  if (!missing(higher.rank)) {
    nodelabels(names(mp.hr), mp.hr, frame = "none", adj = c(1, 1.1), cex = .75)
  }
  # nodelabels(); tiplabels()
  #bubble(mm, phy) # representation by markers
  #append2tips(phy, align = phyUE, grid = TRUE, pch = 19)
  
  if (!missing(min.support) & !is.null(phy$node.label)) 
    node.support(phy$node.label, cutoff = min.support, cex = .6)
  add.scale.bar(length = .05, lwd = 1.5)
  if ( !missing(title) ) title(main = title)
  # legend("topleft", inset = c(.2, .35), legend = names(tax), 
  #        fill = clrs, ncol = 3)
  dev.off()
  invisible(lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv))
}



