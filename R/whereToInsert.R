## This code is part of the megaptera package
## © C. Heibl 2014 (last update: 2018-07-13)

#' @export

whereToInsert <- function(phy, tax, tip, ignore.monophyly, quiet = FALSE){
  
  tip.rank <- "species"
  
  ## Identify anchor point of lowest possible rank
  ## ---------------------------------------------
  lin <- taxdumpLineage(tax, tip)
  if (is.null(lin)) return(NULL)
  for (i in 2:nrow(lin)){
    clade <- lin$taxon[i]
    des <- taxdumpChildren(tax, taxon = clade, tip.rank = tip.rank)
    des <- des$taxon[des$rank == tip.rank]
    des <- setdiff(des, tip)
    des <- intersect(des, phy$tip.label)
    if (length(des)) break
  }
  if (!length(des)) {
    warning("no anchor point for '", tip, "'")
    return(phy)
  }
  an <- noi(phy, des)
  
  
  ## Check monophyly (if clade has more than one tip)
  ## ------------------------------------------------
  lms <- FALSE
  monophyletic <- TRUE
  oops <- setdiff(descendants(phy, an, labels = TRUE), des)
  if (length(oops) & length(des) > 1){
    monophyletic <- FALSE
    ## Find the largest monophyletic subset ('lms') of 'des'
    ## This could be moved upstream to 'noi'
    ## -------------------------------------
    if (!ignore.monophyly){
      test <- extract.clade(phy, an)
      n <- Ntip(test) + (2:Nnode(test))
      lms <- lapply(n, descendants, phy = test, label = TRUE)
      lms <- lms[sapply(lms, function(a, b) all(a %in% b), b = des)]
      ## Empty lms means lms consist of only one species; in this case
      ## we keep the original anchor node 'an'
      ## -------------------------------------
      if (length(lms)){
        lms <- lms[which.max(sapply(lms, length))]
        if (length(lms) > 1) stop("implement me!")
        lms <- unlist(lms)
        an <- noi(phy, lms)
        lms <- TRUE
      } 
    }
  } 
  
  ## Screen output
  ## -------------
  if (!quiet) {
    if (is.null(an)){
      cat("\n", tip, "not in taxonomy")
    } else {
      monophyletic <- ifelse(monophyletic, "(monophyletic)", "(not monophyletic)")
      lms <- ifelse(lms, " (LMS)", "")
      cat("\n", tip, " > anchor clade: ", clade, " ", monophyletic, " > anchor node: ", an, lms, sep = "")
    }
  }
  an
}


