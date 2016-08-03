## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-07-28)

dbReadTaxonomy <- function(megapteraProj, subset, tag){
  
  conn <- dbconnect(megapteraProj)
  
  if ( !dbExistsTable(conn, "taxonomy") )
    stop("no taxonomy table - see ?dbUpdateTaxonomy for help")
  
  ## read taxonomy table
  if ( missing(tag) ){
    tax <- dbGetQuery(conn, "SELECT * FROM taxonomy")
  } else {
    SQL <- paste("SELECT * INTO tmp",
                 "FROM taxonomy",
                 "WHERE", wrapSQL(tag, "tag"))
    dbSendQuery(conn, SQL)
    dbSendQuery(conn, "ALTER TABLE tmp DROP COLUMN tag")
    tax <- dbGetQuery(conn, "SELECT * FROM tmp")
    dbSendQuery(conn, "DROP TABLE tmp")
  }
  
  ## subsetting taxonomy ..
  ## ----------------------
  if ( !missing(subset) ){
    ## .. based on sequence names or ..
    ## --------------------------------
    if ( inherits(subset, "DNAbin") ){
      if ( is.list(subset) ) sset <- names(subset)
      if ( is.matrix(subset) ) sset <- rownames(subset)
      tax <- tax[tax[, megapteraProj@taxon@tip.rank] %in% sset, ]
    }
    ## .. based on sequence names or ..
    ## --------------------------------
    if ( inherits(subset, "phylo") ){
      tax <- tax[tax[, megapteraProj@taxon@tip.rank] %in% subset$tip.label, ]
    }
    ## .. based on <spec.*>
    ## --------------------
    if ( is.character(subset) & length(subset) == 1 ){
      tip.rank <- gsub("_.+$", "", subset)
      subset <- paste("SELECT", tip.rank, 
                      "FROM", subset, 
                      "WHERE status !~ 'excluded'")
      subset <- dbGetQuery(conn, subset)[, tip.rank]
      tax <- tax[tax[, tip.rank] %in% subset, ]
    }
    ## .. based on a vector of species names
    ## -------------------------------------
    if ( is.character(subset) & length(subset) > 1 ){
      tax <- tax[tax[, megapteraProj@taxon@tip.rank] %in% subset, ]
    }
    
    
    ## .. delete ranks lower than 'tip.rank'
    ## -------------------------------------
    if ( megapteraProj@taxon@tip.rank != "spec" ){
      tax <- unique(tax[, 1:which(names(tax) == megapteraProj@taxon@tip.rank)])
    }
  }
  dbDisconnect(conn)
  
  ## remove trailing white spaces
  ## happens often when people prepare taxon lists in Excel
  ## ------------------------------------------------------
  #   tws <- grep(" $|_$", tax[, "spec"])
  #   if ( length(tws) > 0 ){
  #     tax[, "spec"] <- gsub(" $|_$", "", tax[, "spec"])
  #     warning("trailing white space removed in", 
  #             paste("\n - ", head(tax[, "spec"][tws], 3), sep = ""), 
  #             "\n - and ", length(tws) - 3, " more species")    
  #   }
  
  tax <- sqlTaxonomyHeader(tax) # should be obsolete
  tax
}