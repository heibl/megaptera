giConsensus <- function(x, tag = ".", latex = FALSE){
  
  conn <- dbconnect(x)
  
  ## get names of acc_<locus> table
  ## ------------------------------
  acc <- paste("SELECT table_name",
                    "FROM information_schema.tables",
                    "WHERE table_schema='public'",
                    "AND table_type='BASE TABLE'")
  acc <- dbGetQuery(conn, acc)$table_name
  acc <- acc[grep("^acc_", acc)]
  
  
  
  cat("found ", length(acc), " database tables:\n", 
      paste(acc, collapse = ", "), sep = "")
  
  TAG <- paste("AND", wrapSQL(tag, term = "genom"))
  
  obj <- data.frame()
  for ( i in acc ){
    ## check for occureence of 'tag'
    gi <- paste("SELECT DISTINCT genom",
                "FROM", i)
    gi <- dbGetQuery(conn, gi)
    if ( length(grep(tag, gi)) == 0 ) next
    gi <- paste("SELECT taxon, gi",
                "FROM", i,
                "WHERE status ~ 'single|aligned'",
                TAG,
                "ORDER BY taxon, gi")
    gi <- dbGetQuery(conn, gi)
    gi <- data.frame(locus = gsub("^acc_", " ", i), gi,
                     stringsAsFactors = FALSE)
    obj <- rbind(obj, gi)
  }
  dbDisconnect(conn)
  
  stats <- tapply(obj$gi, list(obj$locus, obj$taxon), length)
  stats <- as.vector(stats)
  stats <- stats[!is.na(stats)]
  cat("\ntotal number of sequences:", nrow(obj))
  cat("\nnumber of single sequences:", length(stats[stats == 1]))
  cat("\nnumber of sequences used for conspecific consensus sequences:", 
      sum(stats[stats > 1]))
  cat("\nmedian number of sequences per consensus sequence:", 
      median(stats[stats > 1]))
  
  if ( !missing(latex) ){
    
    obj <- split(obj[, c("taxon", "gi")], obj$locus)
    ff <- function(z){
      z <- split(z$gi, z$taxon)
      z <- sapply(z, paste, collapse = ", ")
      
      spec <- do.call(rbind, strsplit(names(z), "_"))
      d <- duplicated(spec[, 1])
      spec[d, 1] <- paste(substr(spec[d, 1], 1, 1), ".", sep = "")
      spec <- paste(spec[, 1], spec[, 2], sep = "~")
      spec <- paste("\\textit{", spec, "}:", sep = "")
      z <- paste(spec, z)
      paste(z, collapse = "; ")
    }
    obj <- sapply(obj, ff)
    names(obj) <- paste("\\paragraph{", names(obj), "}", sep = "")
    obj <- paste(names(obj), " ", obj, ".", sep = "")
    write(obj, latex)
    
  } else {
    return(obj)
  }
}

