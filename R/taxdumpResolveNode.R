
taxdumpResolveNode <- function(x, y, taxon){
  
  cat(silver(taxon %+% " ... "))
  
  ## root of 'y' will not be considered
  if (taxdump_isRoot(y, taxon)){
    return(x)
  }
  nothing_to_do <- rep(TRUE, 3)
  
  yy <- y[y$taxon == taxon, ]
  if (nrow(yy) > 1) stop("debug me!")
  xx <- x[x$taxon == yy$taxon & x$rank == yy$rank, ]
  if (nrow(xx) > 1) stop("debug me!")
  nothing_to_do[1] <- ifelse(!nrow(xx), FALSE, TRUE) 
  if (nothing_to_do[1]){
    nothing_to_do[2] <- ifelse(xx$status != yy$status, FALSE, TRUE) 
  }
 
  
  if (nothing_to_do){
    cat(green("OK\n"))
    return(x)
  } else {
    repeat{
      yyy <- unique(rbind(yy, y[y$id %in% yy$id, ]))
      xxx <- unique(rbind(xx, x[x$id %in% xx$id, ]))
      if (nrow(yyy) > nrow(yy) | nrow(xxx) > nrow(xx)){
        stop("implement me")
      } else {
        break
      }
    }
    stop("implement me!")
  }
  
  
  
}