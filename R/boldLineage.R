## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-04-11)

#' @importFrom bold bold_tax_id bold_tax_name
#' @export

boldLineage <- function(taxon){
  
  # cat("\n", taxon)
  underscore <- length(grep("_", taxon))
  if (underscore) taxon <- gsub("_", " ", taxon)
  
  obj <- bold_tax_name(taxon)$taxid
  if (is.null(obj)) return(obj)
  obj <- bold_tax_id(obj, includeTree = TRUE)
  obj <- obj[, c("parentid", "taxid", "taxon", "tax_rank")]
  names(obj) <- c("parent_id", "id", "taxon", "rank")
  
  ## order from lower to higher rank
  ## (BOLD does not always return right order)
  ## -----------------------------------------
  id <- which(obj$taxon == taxon)
  while (length(id) < nrow(obj)){
    id <- c(which(obj$id == obj$parent_id[id[1]]), id)
  }
  if (underscore) obj$taxon <- gsub(" ", "_", obj$taxon)
  obj[rev(id), ]
}