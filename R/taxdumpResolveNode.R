#' @export

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
 
  if (all(nothing_to_do)){
    cat(green("OK\n"))
    return(x)
  } else {
    
    ## A: Node is missing
    ## ------------------
    if (!nothing_to_do[1]){
      
      ## Expand nodes
      ## ------------
      repeat{
        yyy <- unique(rbind(yy, y[y$id %in% yy$id, ]))
        xxx <- unique(rbind(xx, x[x$id %in% xx$id, ]))
        if (nrow(yyy) > nrow(yy) | nrow(xxx) > nrow(xx)){
          stop("implement me")
        } else {
          break
        }
      }
      
      ## Double check if taxon is present (might be at another rank, etc ..)
      ## -------------------------------------------------------------------
      if (taxon %in% x$taxon) {
        stop("debug me!")
      }
      
      ## Find anchor point and add node to 'x'
      ## -------------------------------------
      yyLIN <- taxdumpLineage(y, yy$taxon)
      xxLINyy <- intersect(yyLIN$taxon, x$taxon)
      xxLINyy <- xxLINyy[xxLINyy != "root"]
      anchor <- xxLINyy[which.min(match(xxLINyy, yyLIN$taxon))]
      x <- taxdumpAddNode(x, rank = yy$rank, taxon = yy$taxon, parent = anchor)
      cat(red("inserted at node '" %+% bold(anchor) %+% "'\n"))
    }
    
    ## B: Node is of different status
    ## ------------------------------
    if (!nothing_to_do[2]){
      stop("implement me [B]")
    }
    
    ## C: Lineage is different
    ## -----------------------
    if (!nothing_to_do[3]){
      stop("implement me [C]")
    }

  }
  
  
  x
}