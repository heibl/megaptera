## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-03-24)

#' @export

taxdumpDaughters <- function(x, taxon, tip.rank, indet = indet.strings()){
  
  ## checks
  ## ------
  if (length(taxon) > 1) stop("cannot handle more than one taxon")
  
  ## get taxonomy if necessary (i.e. x is not a parent-child-table)
  ## --------------------------------------------------------------
  if (inherits(x, "megapteraProj")){
    x <- dbReadTaxonomy(x)
  }
  
  ## get id of 'taxon', i.e. the root node
  ## -------------------------------------
  if (length(grep("[[:alpha:]]", taxon)) == 1){
    id <- x[x$taxon == taxon, "id"]
    if (!length(id)) {
      warning("taxon '", taxon, "' not available")
      return(NULL)
    }
  } else {
    id <- taxon
  }
  
  ## get all daughter nodes
  ## ----------------------
  this.id <- id
  gain <- length(id)
  while (gain > 0){
    this.id <- x[x$parent_id %in% this.id, "id"]
    id <- c(id, this.id)
    gain <- length(this.id)
  }
  x <- x[x$id %in% id, c(2, 1, 4, 3)]
  rownames(x) <- NULL
  
  ## remove *species* that do not correspond to structurally
  ## valid Latin binomials according to 'indet'
  ## ------------------------------------------
  notvalid <- grep(indet.strings(collapse = TRUE), x$taxon)
  notvalid <- intersect(notvalid, which(x$rank == "species"))
  if (length(notvalid) > 0){
    # message(length(notvalid), " taxa removed")
    x <- x[-notvalid, ]
  }
  
  ## remove rows of rank below 'tip.rank'
  ## ------------------------------------
  id <- vector()
  this.id <- x[x$rank == tip.rank, "id"]
  gain <- length(this.id)
  while (gain > 0){
    this.id <- x[x$parent_id %in% this.id, "id"]
    id <- c(id, this.id)
    gain <- length(this.id)
  }
  x <- x[!x$id %in% id, ]
  
  x
}