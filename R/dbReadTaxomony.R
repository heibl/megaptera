## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-11-09)

#' @rdname dbTaxonomy
#' @import DBI
#' @export

dbReadTaxonomy <- function(megProj, tip.rank, subset, tag, root = "tol"){
  
  if (missing(tip.rank)) tip.rank <- megProj@taxon@tip.rank
  root <- match.arg(root, c("tol", "mrca"))
  
  conn <- dbconnect(megProj)
  
  if (!dbExistsTable(conn, "taxonomy"))
    stop("no taxonomy table - see ?dbUpdateTaxonomy for help")
  
  ## read taxonomy table
  ## -------------------
  if (missing(tag)){
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
  
  ## truncate taxonomy to tip rank
  ## -----------------------------
  id <- tax[tax$rank == tip.rank, "id"]
  tdDescendants <- function(tax, id){
    all_ids <- vector()
    gain <- length(id)
    while (gain > 0){
      id <- tax[tax$parent_id %in% id, "id"]
      all_ids <- c(all_ids, id)
      gain <- length(id)
    }
    all_ids
  }
  id <- lapply(id, tdDescendants, tax = tax)
  id <- unlist(id)
  tax <- tax[!tax$id %in% id, ]
  
  ## subsetting taxonomy ..
  ## ----------------------
  if (!missing(subset)){
    
    ## .. based on sequence names or ..
    ## --------------------------------
    if (is.character(subset)) {
      ## .. based on MAS table or ..
      ## ---------------------------
      if (length(grep("^[genus]|[species]_", subset))){
        tip.rank <- gsub("_.+$", "", subset)
        subset <- paste("SELECT taxon", 
                        "FROM", subset, 
                        "WHERE status !~ 'excluded'")
        subset <- dbGetQuery(conn, subset)[, "taxon"]
      } else {
        ## .. based on a vector of species names (without synonyms)
        ## --------------------------------------------------------
        subset <- subset
      }
    } 
    
    ## .. based on a phylogeny or ..
    ## -----------------------------
    if (inherits(subset, "phylo")){
      subset <-  subset$tip.label
    } else {
      if (inherits(subset, "DNAbin")){
        if (is.list(subset)) subset <- names(subset)
        if (is.matrix(subset)) subset <- rownames(subset)
      } else {
        ## .. based on a list of species names (with synonyms)
        ## ---------------------------------------------------
        if (is.list(subset)) subset <- sapply(subset, head, n = 1)
      }
    }
    
    ## tip rank is genus, but subset species-level
    ## -------------------------------------------
    if (tip.rank == "genus" & all(is.Linnean(subset))){
      subset <- unique(strip.spec(subset))
    }
    
    ## do the actual subsetting
    ## ------------------------
    subset <- gsub("_", " ", subset)
    id <- all_ids <- tax[tax$taxon %in% subset, "id"]
    while (length(id) > 1) {
      id <- unique(tax[tax$id %in% id, "parent_id"])
      all_ids <- c(all_ids, id)
    }
    tax <- tax[tax$id %in% all_ids, ]
    
    ## test if taxa that are present in subset
    ## are missing from taxonomy table
    ## -------------------------------
    id <- setdiff(subset, tax$taxon)
    if (length(id)){
      warning(length(id), " taxa of 'subset' are missing from taxonomy table:",
              paste("\n", sort(id), collapse = ""))
    }
  }
  
  ## remove tree-of-life-root (if necessary, see 2nd condition)
  ## ----------------------------------------------------------
  if (root == "mrca" & 1 %in% tax$parent_id){
    nodes <- 1
    repeat{
      nn <- tax[tax$parent_id == nodes[1], "id"]
      nn <- nn[nn != 1]
      if (length(nn) > 1) break
      nodes <- c(nn, nodes)
    }
    tax <- tax[!tax$id %in% nodes[-1], ]
  }
  dbDisconnect(conn)
  tax
}