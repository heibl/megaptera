## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-02-22)

megaptera2Rmarkdown <- function(x, file){
  
  if ( missing(file) ) file <- "megaptera.Rmd"
  
  z <- c("---", 
         paste('title: "', 
               paste("megaptera", packageDescription("megaptera")$Version),
               '"', sep = ""),
         paste("date:", Sys.time()),
         "output:", 
         "  html_document:",
         "    toc: true",
         paste("    css:", gsub("Rmd", "css", file)),
         "---")
  
  save(x, file = "megapteraProj.rda")
  
  z <- c(z, 
         "```{r, echo=FALSE, message=FALSE, warning=FALSE}",
         "library(megaptera)",
         "load('megapteraProj.rda')", 
         "```")
  
  ## SETTINGS: TAXONOMY
  ## ------------------
  spec.list <- all(unlist(sapply(x@taxon@ingroup, is.Linnean)))
  ingroup <- head(x@taxon@ingroup, 3)
  if ( spec.list ) ingroup <- paste("*", ingroup, "*", sep = "")
  ingroup <- paste(ingroup, collapse = ", ")
  n <- length(x@taxon@ingroup)
  if ( n > 3 )
    ingroup <- paste(ingroup, paste(", ... [", n, " taxa in total]", 
                                    sep = ""), sep = "")
  outgroup <- head(x@taxon@outgroup, 3)
  outgroup <- paste("*", outgroup, "*", sep = "")
  outgroup <- paste(outgroup, collapse = ", ")
  n <- length(x@taxon@outgroup)
  if ( n > 3 )
    outgroup <- paste(outgroup, paste(", ... [", n, " taxa in total]", 
                                      sep = ""), sep = "")
  z <- c(z, 
         "# Settings",
         "## Taxonomic settings",
         paste("- **Kingdom**:", x@taxon@kingdom),
         paste("- **Ingroup**:", ingroup),
         paste("- **extend ingroup**:", ifelse(x@taxon@extend.ingroup, "yes", "no")),
         paste("- **Outgroup**:", outgroup),
         paste("- **extend outgroup**:", ifelse(x@taxon@extend.outgroup, "yes", "no")),
         paste("- **Hybrids**:", ifelse(x@taxon@hybrids, "included", "excluded")),
         paste("- **Tip rank**:", x@taxon@tip.rank),
         paste("- **Reference rank** :", x@taxon@reference.rank),
         paste("- **User-defined guide tree** :", ifelse(inherits(x@taxon, "taxonGuidetree"), "yes", "no")),
         "")
  
  ## SETTINGS: DATABASE
  ## ------------------
  z <- c(z, 
         "## PostgreSQL connection parameters",
         paste("- **host**:", x@db@host),
         paste("- **port**:", x@db@port),
         paste("- **dbname**:", x@db@dbname),
         paste("- **user**:", x@db@user),
         paste("- **password**:", x@db@password),
         "")
  
  ## SETTINGS: PIPELINE
  ## ------------------
  z <- c(z, 
         "## Pipeline settings",
         paste("- **Execution**:", ifelse(x@params@parallel, "parallel", "serial")),
         paste("- **Number of CPUs**:", x@params@cpus),
         paste("- **Type of cluster**:", x@params@cluster.type),
         "")
  
  ## 2: STATUS locus-wise
  ## -----------------
  tabs <- dbTableNames(x)
  tabs <- tabs[grep("acc_", tabs)]
  if ( length(tabs) == 0 ){
    ss <- c(TRUE, rep(FALSE, 8))
    names(ss) <- LETTERS[1:9]
    ss <- do.call(rbind, list(ss))
    ss <- htmlTable(ss)
    z <- c(z, "# Status of the pipeline", ss, "-")
    write(z, file = file)
    return()  
  } 
  tabs <- gsub("acc_", "", tabs)
  STATUS <- lapply(tabs, checkStatus, x = x)
  STATUS <- do.call(rbind, STATUS)
  rownames(STATUS) <- tabs
  id <- order(apply(STATUS, 1, which.min),
              decreasing = TRUE)
  STATUS <- STATUS[id, , drop = FALSE]
  
  ss <- cbind(rownames(STATUS), STATUS)
  ss[ss == TRUE] <- '<div class="gd"></div>'
  ss[ss == "FALSE"] <- '<div class="rd"></div>'
  ss <- htmlTable(ss)
  z <- c(z, 
         "# Status of the pipeline",
         ss, "")
  
  ## 3: TAXONOMY
  ## ------------------
  tax <- dbReadTaxonomy(x)
  taxonomy2html(tax)
  tax$spec <- gsub("_", " ", tax$spec)
  ingroup.queried <- x@taxon@ingroup
  if ( spec.list ){
    ingroup.queried <- sapply(ingroup.queried, head, n = 1)
    ingroup.received  <- tax$spec[grep("^ingroup", tax$tag)]
    extended.ingroup.received <- tax$spec[grep("extended ingroup", tax$tag)]
    ingroup.missing <- setdiff(ingroup.queried, ingroup.received)
    ingroup.false <- setdiff(ingroup.received, ingroup.queried)
    if ( length(ingroup.false) > 0 ){
      ingroup.false <- paste("- **WARNING**: *", ingroup.false, 
                             "* falsely as ingroup classified", sep = "")
    } else {
      ingroup.false <- ""
    }
    tax.coverage <- c(
      "Query of the NCBI Taxonomy Database gave these results:", "",
      paste("- **", length(ingroup.received), "** out of ", 
            length(ingroup.queried), " ingroup species (", 
            round(length(ingroup.received) / length(ingroup.queried) * 100, 2),
            "%) were found", sep = ""),
      paste("- **", length(ingroup.missing), "** ingroup species (", 
            round(length(ingroup.missing) / length(ingroup.queried) * 100, 2),
            "%) were missing", sep = ""),
      paste("- **", length(extended.ingroup.received), 
            "** congenerics of the ingroup species were included in the extended ingroup", 
            sep = ""),
      ingroup.false)
    
  } else {
    tax.coverage <- c(
      "Query of the NCBI Taxonomy Database gave these results:", "",
      paste("- **", length(grep("ingroup", tax$tag)), 
            "** ingroup species", sep = ""),
      paste("- **", length(grep("outgroup", tax$tag)), 
            "** outgroup species", sep = "")
    )
  }
  ## guide tree
  ## ----------
  z <- c(z, "# Taxonomy",
         tax.coverage, "",
         "<a href='taxonomy.html'>view taxonomy table</a>",
         "",
         #          "```{r, echo=FALSE, message=FALSE}",
         #          "if ( inherits(x@taxon, 'taxonGuidetree') ){",
         #          "gt <- comprehensiveGuidetree(x, tip.rank = 'gen')",
         #          "} else {" ,
         #          "gt <- dbReadTaxonomy(x)",
         #          "gt <- tax2tree(gt, tip.rank = 'gen')",
         #          "}",
         #          "gt <- ladderize(gt)",
         #          "plot(gt, type = 'clado')",
         # "```", 
         "")
  
  ## RAW SEQUENCES
  ## -------------
  z <- c(z, "# DNA sequence retrieval",
  "## Species without any sequences",
  "## Number of species per locus",
  "```{r, echo=FALSE, message=FALSE}",
  "gb <- checkSpecLocus(x, 'gb')",
  "```", "")


## SELECTED SEQUENCES
## ------------------
if ( all(!STATUS[, "F"]) ){
  write(z, file = file)
  return()  
}
tab <- checkSpecLocus(x, 'sel', plot = FALSE)
tab <- cbind(rownames(tab$specPerMarker), tab$specPerMarker)
colnames(tab) <- c("Locus", "Number of sequences", 
                   "Number of private sequences")
tab <- htmlTable(tab) 
z <- c(z, 
       "# Sequence selection and alignment",
       "## Species not selected for alignment",
       "## Number of species per locus",
       "```{r, echo=FALSE, message=FALSE}",
       "sel <- checkSpecLocus(x, 'sel')",
       "```", "", tab)

## SATURATION + BLOCKS
## ------------------
msa <- dbSummaryMSA(x)
msa[, -1] <- apply(msa[, -1], 2, 
                   function(y) round(as.numeric(y), digits = 4))
msa <- cbind(Locus = rownames(msa), msa)
tab <- htmlTable(msa)
z <- c(z, 
       "# Evaluation of alignments",
       "## Saturation",
       tab, 
       "```{r, echo=FALSE, message=FALSE}",
       "b <- checkBlocks(x)",
       "```", "")




write(z, file = file)
}