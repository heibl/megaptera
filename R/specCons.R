## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-03-28)

#' @export

specCons <- function(obj, log) {
  
  if (!missing(log)){
    acc <- splitGiTaxon(rownames(obj))
    acc <- split(acc$gi, acc$taxon)
    # gen <- ifelse(length(acc) > 1), strip.spec(names(acc[1])), NULL)
    acc <- sapply(acc, paste, collapse = ", ")
    acc <- c(paste("\n Genus:", strip.spec(names(acc[1]))),
             paste("\n >", names(acc), "-", acc))
    slog(acc, file = log)
  }
  
  if (nrow(obj) > 1){
    obsvalue <- as.raw(c(136, 40, 72, 24))
    tableRaw <- function(obj, obsvalue){
      sapply(obsvalue, function(o, obj) length(which(obj %in% o)), obj = obj)
    }
    ## corresponds to method 'profile' in seqinr::consensus
    obj <- apply(obj, 2, tableRaw, obsvalue = obsvalue)
    return(apply(obj, 2, function(obj, o) o[which.max(obj)], o = obsvalue))
  } else {
    return(as.vector(obj))
  }
}