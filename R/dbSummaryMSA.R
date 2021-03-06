## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2018-01-31)

#' @importFrom stats mad sd
#' @export

dbSummaryMSA <- function(x, masked = FALSE){
  
  ## Which loci completed stepH successfully?
  ## ----------------------------------------
  conn <- dbconnect(x)
  loci <- dbGetQuery(conn, "SELECT locus, step_h FROM progress")
  loci <- loci$locus[loci$step_h == "success"]
  dbDisconnect(conn)
  
  ## Core Function
  ## -------------
  engine <- function(x, locus, model = "JC69"){
    
    s <- dbReadMSA(x, locus)
    if ( is.null(s) ) return(NULL)
    p.uncorr <- dist.dna(s, model = "raw", pairwise.deletion = TRUE)
    p.corr <- dist.dna(s, model = model, pairwise.deletion = TRUE)
    obj <- list(
      n = nrow(s),
      min = min(p.uncorr, na.rm = TRUE),
      max = max(p.uncorr, na.rm = TRUE),
      median = median(p.uncorr, na.rm = TRUE),
      sd = sd(p.uncorr, na.rm = TRUE),
      mad = mad(p.corr - p.uncorr, na.rm = TRUE)
    )
    obj
  }
  
  ## Assemble information
  ## --------------------
  obj <- lapply(loci, engine, x = x)
  names(obj) <- loci
  obj <- obj[!sapply(obj, is.null)]
  obj <- obj[order(sapply(obj, function(z) z$mad))]
  obj <- do.call(rbind, obj)
  rownames(obj) <- gsub(paste(x@taxon@tip.rank, "_", sep = ""), "", rownames(obj))
  rownames(obj) <- gsub("_", "", rownames(obj))
  as.matrix(obj)
}