dbSummaryMSA <- function(x, masked = FALSE){
  
  tabs <- dbTableNames(x, x@taxon@tip.rank)
  
  engine <- function(msa, masked, model = "JC69"){
    
    s <- dbReadDNA(x, msa, taxon = ".*", regex = TRUE,
                   ignore.excluded = TRUE, masked = masked, 
                   blocks = "ignore")
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
  obj <- lapply(tabs, engine, masked = masked)
  names(obj) <- tabs
  obj <- obj[!sapply(obj, is.null)]
  obj <- obj[order(sapply(obj, function(z) z$mad))]
  obj <- do.call(rbind, obj)
  rownames(obj) <- gsub(paste(x@taxon@tip.rank, "_", sep = ""), "", rownames(obj))
  rownames(obj) <- gsub("_", ".", rownames(obj))
  as.matrix(obj)
}