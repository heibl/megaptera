## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-07-05)

stepA <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  
  ## iniate logfile
  ## --------------
  logfile <- "stepA.log"
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""), 
       "\nSTEP A: searching and downloading taxonomy from GenBank\n",
       file = logfile)
  
  if ( file.exists("ncbiTaxonomy-missing.txt") ) unlink("ncbiTaxonomy-missing.txt")
  
  conn <- dbconnect(x)
  notable <- !dbExistsTable(conn, "taxonomy")
  dbDisconnect(conn)
  
  if (  notable | x@update ){
    
    ## download NCBI taxonomy for ingroup
    ingroup <- ncbiTaxonomy(taxon = x@taxon@ingroup, 
                            extend = x@taxon@extend.ingroup,
                            kingdom = x@taxon@kingdom, 
                            megapteraProj = x)
    if ( is.null(ingroup) ){
      stop("could not retrieve taxonomy for ingroup")
    }
    ingroup <- fixTaxonomy(ingroup, auto = TRUE, ignore = "synonym")
    tag <- rep("ingroup (NCBI)", nrow(ingroup))
    if ( x@taxon@extend.ingroup ){
      extended.species.list <- !ingroup$spec %in% sapply(x@taxon@ingroup, head, 1)
      if ( any(extended.species.list) ){
        tag[extended.species.list] <- "extended ingroup (NCBI)"
      }
    }
    dbUpdateTaxonomy(x, ingroup, tag = tag)
    
    ## download NCBI taxonomy for ingroup
    outgroup <- ncbiTaxonomy(taxon = x@taxon@outgroup, 
                             kingdom = x@taxon@kingdom, 
                             megapteraProj = x)
    if ( is.null(outgroup) ){
      stop("could not retrieve taxonomy for outgroup")
    }
    outgroup <- fixTaxonomy(outgroup, auto = TRUE, ignore = "synonym")
    dbUpdateTaxonomy(x, outgroup, tag = "outgroup (NCBI)")
    
  } else {
    slog("\ntaxonomy already downloaded", file = logfile)
  }
  
  slog("\n\nSTEP A finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
}