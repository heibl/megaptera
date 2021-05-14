## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2021-03-18)

#' @rdname dbTaxonomy
#' @import DBI
#' @export

dbReadTaxonomy <- function(megProj, tip.rank, subset, tag, root = "tol", 
                           syn = FALSE, drop.extinct = TRUE){
  
  ## CHECKs
  if (!inherits(megProj, "megapteraProj"))
    stop("'megProj' is not of class 'megapteraProj'")
  if (!missing(subset)){
    if (megProj@locus@kind == "undefined" & length(subset) == 1) ## i.e. subset == "species_sequence"
      stop("locus must be defined for subsetting based on table 'species_sequence'")
  }
  
  if (missing(tip.rank)) tip.rank <- megProj@taxon@tip.rank
  root <- match.arg(root, c("tol", "mrca"))
  
  conn <- dbconnect(megProj)
  
  if (!dbExistsTable(conn, "taxonomy"))
    stop("no taxonomy table - see ?dbUpdateTaxonomy for help")
  
  ## Read taxonomy table
  ## -------------------
  if (missing(tag)){
    SQL <- "SELECT * FROM taxonomy"
    if (!syn) SQL <- paste(SQL, "WHERE status != 'synonym'")
    tax <- dbGetQuery(conn, SQL)
  } else {
    if (syn) stop("implement me!")
    SQL <- paste("SELECT * INTO tmp",
                 "FROM taxonomy",
                 "WHERE", wrapSQL(tag, "tag"))
    dbSendQuery(conn, SQL)
    dbSendQuery(conn, "ALTER TABLE tmp DROP COLUMN tag")
    tax <- dbGetQuery(conn, "SELECT * FROM tmp")
    dbSendQuery(conn, "DROP TABLE tmp")
  }
  
  ## Truncate taxonomy to tip rank; this is done in two steps:
  ## ---------------------------------------------------------
  
  ## Step 1: Remove all nodes below tip.rank
  ## ---------------------------------------
  id <- tax$id[tax$rank == tip.rank & tax$status == "scientific name"]
  tdDescendants <- function(tax, id){
    all_ids <- vector()
    gain <- length(id)
    while (gain > 0){
      id <- tax$id[tax$parent_id %in% id & tax$status == "scientific name"]
      all_ids <- c(all_ids, id)
      gain <- length(id)
    }
    all_ids
  }
  id <- lapply(id, tdDescendants, tax = tax)
  id <- unlist(id)
  tax <- tax[!tax$id %in% id, ]
  
  ## Step 2
  ## There can be lineages with tip.rank missing,
  ## e.g. subgenus: Neocicindela, no genus, tribe: Cicindelini,
  ## These lineages will be dropped entirely.
  ## ----------------------------------------
  if (drop.extinct){
    tn <- taxdump_isTerminal(tax)
    id <- tax$id[tn & tax$rank != tip.rank & tax$status != "synonym"]
    if (length(id)){
      if (megProj@params@debug.level >= 3){
        warning(length(id)," terminal taxa without a taxon of rank '", tip.rank, 
                "' in their lineage were removed:",
                paste("\n-", tax$taxon[tax$id %in% id]))
      }
      tax <- taxdumpDropTip(tax, id)
    }
  }
  
  
  ## Subsetting taxonomy ..
  ## ----------------------
  if (!missing(subset)){
    
    ## .. based on sequence names or ..
    ## --------------------------------
    if (is.character(subset)) {
      ## .. based on MSA table or ..
      ## ---------------------------
      if (length(grep("^[genus]|[species]_", subset)) == 1){
        tip.rank <- gsub("_.+$", "", subset)
        subset <- paste("SELECT taxon", 
                        "FROM", subset, 
                        "WHERE status !~ 'excluded'",
                        "AND", wrapSQL(megProj@locus@sql, "locus", "="))
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
  
  ## Remove tree-of-life-root (if necessary, see 2nd condition)
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