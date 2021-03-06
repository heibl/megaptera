## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2018-06-26)

#' @title Add Tips to a Phylogenetic Tree
#' @description Add tips (species) to a phylogenetic tree according to their
#'   taxonomic classification.
#' @param phy An object of class \code{\link{phylo}}.
#' @param tip A character string giving the name of the species to be added.
#' @param sister A character string giving the name of the species or genus to
#'   which \code{tip} shuold be added as sister; this argument overrides
#'   \code{tax}.
#' @param tax A data frame containing the taxonomic classification.
#' @param insert A character string indicating the positions where the species
#'   is to be inserted: \code{"crown"}, \code{"stem"}, \code{"randomly"}, or any
#'   unambigous abbreviation of these. This option will only have an effect if
#'   \code{phy} contains more than one congeneric of \code{tip}.
#' @param ignore.monophyly Logical, indicating if monophyly should be considered
#'   when selecting the anchor point for \code{tips}.
#' @param stem.edge A real number greater than 0 and smaller than 1, which gives
#'   the fraction of the terminal branch length that will be assigned to the
#'   branch subtending the newly created MRCA of \code{tip} and its single
#'   congeneric species in \code{phy}. Will have no effect if \code{phy}
#'   contains more than one congeneric of \code{tip}.
#' @param quiet Logical, indicating if screen output should be suppressed.
#' @return An object of class \code{\link{phylo}}.
#' @seealso \code{\link{dbReadTaxonomy}} for reading a taxonomy table from the
#'   postgreSQL database.
#' @importFrom ape extract.clade Ntip Nnode Nedge rotate is.ultrametric
#'   branching.times
#' @importFrom ips descendants noi fixNodes tipHeights
#' @importFrom utils head
#' @importFrom stats runif
#' @export

 
addSingleTip <- function(phy, tip, sister, tax, 
                         insert = "crown", 
                         ignore.monophyly,
                         stem.edge = 0.5,
                         quiet = FALSE){
  
  ## Check arguments for validity
  ## ----------------------------
  if (!inherits(phy, "phylo")) 
    stop("'phy' is not of class 'phylo'")
  insert <- match.arg(insert, c("crown", "stem", "randomly", "below"))
  
  ## Enforce underscore
  ## ------------------
  tip <- gsub(" ", "_", tip)
  if (missing(sister)){
    tax$taxon[tax$rank == "species"] <- gsub(" ", "_", tax$taxon[tax$rank == "species"])
  }
  
  ## Check if tip is already present in phy
  ## --------------------------------------
  if (tip %in% phy$tip.label) {
    message("'tip' is already contained in 'phy'")
    return(phy)
  }
  
  ## Number of tips, internal nodes and edges:
  ## ----------------------------------------
  nt <- Ntip(phy); ni <- Nnode(phy); ne <- Nedge(phy)
  
  
  ## Identify anchor point
  ## ---------------------
  if (missing(sister)){
    
    ## Case 1: Anchor point is taken from taxonomy
    ## -------------------------------------------
    an <- whereToInsert(phy, tax, tip, ignore.monophyly, quiet)
    if (is.null(an)) {
      if (!quiet) {
        message("no anchor point for '", tip, "'")
      }
      return(phy)
    } 
    
    ## Index 'lower' won't work when tip is inserted
    ## at node nt; rotating subtending node of nt
    ## circumvents this problem
    if (an == nt){
      phy <- fixNodes(rotate(phy, phy$edge[phy$edge[, 2] == nt, 1]))
      an <- whereToInsert(phy, tax, tip)
    }
  } else {
    
    ## Case 2: Anchor point is sister species or sister genus
    ## ------------------------------------------------------
    if (length(grep("_", sister)) == 1){
      an <- which(phy$tip.label == sister)
      if (!length(an)) stop("sister species '", sister, "' not in 'phy'")
      if (insert == "stem"){
        an <- phy$edge[phy$edge[, 2] == an, 1]
        insert <- "crown"
      }
    } else {
      an <- noi(phy, paste0("^", sister, "_"), regex = TRUE)
      insert <- "sister"
    }
  }
  
  if (an <= nt){
    
    ## Add tip to one congeneric
    ## -------------------------
    pretip <- an
    
    ## an: This is now the subtending (or stem) node
    ## of the node where the new tip is to be inserted
    an <- phy$edge[phy$edge[, 2] == pretip, 1]
    
    ## New internal node number
    newinternal <- descendants(phy, an, type = "i")
    if (length(newinternal) == 0){
      newinternal <- an + 2
    } else {
      newinternal <- max(descendants(phy, an, type = "i")) + 2
    }
    
    # Increase node number greater than 'pretip' by 1
    # to create a gap to insert new tip at number pretip + 1
    phy$edge[phy$edge > pretip] <- phy$edge[phy$edge > pretip] + 1
    an <- an + 1 # step ancestral node accordingly
    ## create splits above and below gap
    id <- which(phy$edge[, 1] == an & phy$edge[, 2] == pretip)
    upper <- 1:id; lower <- (id + 1):ne
    phy$edge[phy$edge >= newinternal] <- phy$edge[phy$edge >= newinternal] + 1
    phy$edge[id, 2] <- newinternal
    
    # add edges
    phy$edge <- rbind(phy$edge[upper, ], 
                      c(newinternal, pretip),
                      c(newinternal, pretip + 1),
                      phy$edge[lower, ])
    
    # add edge.lengths
    phy$edge.length <- c(phy$edge.length[head(upper, -1)],
                         phy$edge.length[id] * stem.edge,
                         rep(phy$edge.length[id] * (1 - stem.edge), 2),
                         phy$edge.length[lower])
    
    # add tip.label
    phy$tip.label <- c(phy$tip.label[1:pretip], 
                       tip,
                       phy$tip.label[(pretip + 1):nt])
    
    # ajust internal node number
    phy$Nnode <- ni + 1 # and not nt, which would be valid 
    # only for binary trees!
    
  } else {
    
    if (insert == "sister"){
      
      ## Add tip as sister to another genus (2018-06-26)
      ## -----------------------------------------------
      
      ## 1. Identify relevant edge and tip
      id <- which(phy$edge[, 2] == an)
      new_tip <- max(descendants(phy, an, type = "t"))
      id2 <- which(phy$edge[, 2] == new_tip)
      
      ## 2. Calculate edge length of new tip
      ## -----------------------------------
      new_edge_length <- an
      while (new_edge_length[1] > nt){
        new_edge_length <- c(min(phy$edge[phy$edge[, 1] == new_edge_length[1], 2]), 
                             new_edge_length)
      }
      new_edge_length <- which(phy$edge[, 2] %in% rev(new_edge_length)[-1])
      new_edge_length <- sum(phy$edge.length[new_edge_length])
      
      ## 3. Step node numbers >= an by one
      phy$edge[phy$edge >= an] <- phy$edge[phy$edge >= an] + 1
  
      ## 4. Step node numbers > new_tip by one
      phy$edge[phy$edge > new_tip] <- phy$edge[phy$edge > new_tip] + 1
      
      ## 5. Insert new edges
      upper <- 1:(id - 1)
      sister <- (id + 1):id2
      lower <- (id2 + 1):ne
      phy$edge <- rbind(phy$edge[upper, ], 
                        c(an, an + 1),
                        c(an + 1, an + 2),
                        phy$edge[sister, ],
                        c(an + 1, new_tip + 1),
                        phy$edge[lower, ])
      
      ## 6. Insert new edge lengths
      phy$edge.length <- c(phy$edge.length[upper],
                           phy$edge.length[id] * stem.edge,
                           phy$edge.length[id] * (1 - stem.edge),
                           phy$edge.length[sister],
                           phy$edge.length[id] * (1 - stem.edge) + new_edge_length,
                           phy$edge.length[lower])
      
      ## 7. Insert tip label
      phy$tip.label <- c(phy$tip.label[1:new_tip], 
                         tip,
                         phy$tip.label[(new_tip + 1):nt])
      
      ## 8. Adjust internal node number
      phy$Nnode <- ni + 1
      
    } else {
      
      ## Add tip to more than one congenerics or a higher rank
      ## -----------------------------------------------------    
      if (insert == "stem"){
        # an <- noi(phy, strip.spec(tip), regex = TRUE, stem = TRUE)
        an <- phy$edge[phy$edge[, 2] == an, 1]
      } 
      if (insert == "randomly") {
        an <- descendants(phy, an, type = "i")
        an <- sample(an, 1)
      }
      
      # Crown group age:
      if (is.ultrametric(phy)){
        ael <- branching.times(phy)[an - Ntip(phy)]
      } else {
        tip.heights <- tipHeights(extract.clade(phy, an))
        ael <- runif(1, min(tip.heights), max(tip.heights))
      }
      
      # Number of tip to be inserted: the smallest tip number,
      # because the existing tip number will be increased by 1
      newtip <- min(descendants(phy, an, "t"))
      
      # Increase number of internal nodes by 1
      phy$edge[phy$edge >= newtip] <- phy$edge[phy$edge >= newtip] + 1
      an <- an + 1
      
      # insert before id (after would be more difficult)
      id <- min(which(phy$edge[, 1] == an))
      upper <- 1:(id - 1); lower <- id:ne
      
      # add edge
      phy$edge <- rbind(phy$edge[upper, ], 
                        c(an, newtip),
                        phy$edge[lower, ])
      # add edge.length
      phy$edge.length <- c(phy$edge.length[upper], 
                           ael,
                           phy$edge.length[lower])
      ## ADD tip.label
      ## -------------
      ## three cases need to be treated:
      ## (1) the special case that the new tip will 
      ## be the *first* tip in the tree, 
      ## (2) the special case that the new tip will 
      ## be the *last* tip in the tree,
      ## (3) the general case of the new tip's number 
      ## being intermediate
      if (newtip == 1){
        phy$tip.label <- c(tip,
                           phy$tip.label[newtip:nt])
      } 
      if (newtip == nt){ 
        phy$tip.label <- c(phy$tip.label[1:(newtip - 1)], 
                           tip)
      } 
      if (!newtip %in% c(1, nt)) {
        phy$tip.label <- c(phy$tip.label[1:(newtip - 1)], 
                           tip,
                           phy$tip.label[newtip:nt])
      }
    }
  }
  #   fixNodes(phy)
  phy
}