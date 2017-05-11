## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-05-10)

#' @title RMarkdown Report
#' @description Creates a status report for a megaptera Project using RMarkdown.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param file A character vector giving a file name.
#' @return None, \code{megaptera2Rmarkdown} is called for its side effect.
#' @seealso \code{\link{checkSpecLoc}}, \code{\link{checkSpecies}}
#' @export

megaptera2Rmarkdown <- function(x, file){
  
  if (missing(file)) file <- "megaptera.Rmd"
  if (length(grep("^report/", file)) == 0){
    file <- paste("report", file, sep = "/")
  }
  
  ## Prepare some strings for output
  ## -------------------------------
  tip_rank <- ifelse(x@taxon@tip.rank == "genus", "genera", "species")
  Tip_rank <- ifelse(x@taxon@tip.rank == "genus", "Genera", "Species")
  
  z <- c("---", 
         paste0('title: "', 
                paste("megaptera", packageDescription("megaptera")$Version),
                '"'),
         paste("date:", Sys.time()),
         "output:", 
         "  html_document:",
         "    toc: true",
         "    css: megaptera.css",
         "---")
  
  save(x, file = "report/megapteraProj.rda")
  
  z <- c(z, 
         "```{r, echo=FALSE, message=FALSE, warning=FALSE}",
         "library(megaptera)",
         "load('megapteraProj.rda')", 
         "```")
  
  ## SETTINGS: TAXONOMY
  ## ------------------
  spec.list <- all(unlist(sapply(x@taxon@ingroup, is.Linnean)))
  ingroup <- head(x@taxon@ingroup, 3)
  if (spec.list) ingroup <- paste("*", ingroup, "*", sep = "")
  ingroup <- paste(ingroup, collapse = ", ")
  n <- length(x@taxon@ingroup)
  if (n > 3){
    ingroup <- paste(ingroup, paste(", ... [", n, " taxa in total]", 
                                    sep = ""), sep = "")
  }
  outgroup <- head(x@taxon@outgroup, 3)
  outgroup <- paste("*", outgroup, "*", sep = "")
  outgroup <- paste(outgroup, collapse = ", ")
  n <- length(x@taxon@outgroup)
  if (n > 3){
    outgroup <- paste(outgroup, paste(", ... [", n, " taxa in total]", 
                                      sep = ""), sep = "")
  }
  
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
  tabs <- dbTableNames(x, "acc")
  if (!length(tabs)){
    if ("taxonomy" %in% dbTableNames(x, "all")){
      z <- c(z, "# Status of the pipeline", 
             "Taxonomy table has been created (by running stepA)", "",
             "No sequences have been downloaded yet (stepB has to be called next)", "")
    }
  } else {
    
    ## create status table when stepB has been called
    ## at least once
    STATUS <- dbProgress(x, locus.wise = FALSE)
    names(STATUS) <- c("Locus", LETTERS[2:9])
    ss <- STATUS
    ss[ss == "pending"] <- 1
    ss[ss == "success"] <- 2
    ss[ss == "failure"] <- 3
    ss[ss == "error"] <- 4
    ss <- ss[order(ss[, 2], ss[, 3], ss[, 4], ss[, 5], ss[, 6], ss[, 7], ss[, 8], ss[, 9]), ]
    ss$Locus <- gsub("^_", "", ss$Locus)
    
    ss[ss == 1] <- '<div class="pending"></div>'
    ss[ss == 2] <- '<div class="success"></div>'
    ss[ss == 3] <- '<div class="failure"></div>'
    ss[ss == 4] <- '<div class="error"></div>'
    ss <- htmlTable(ss)
    z <- c(z, 
           "# Status of the pipeline",
           "<table style=\"border-collapse: collapse;\">",
           "<tr>",
           "<td style=\"border:2px solid #fff;\"><div class=\"pending\"></div></td>",
           "<td style=\"border:2px solid #fff;\">pending, </td>",
           "<td style=\"border:2px solid #fff;\"><div class=\"success\"></div></td>",
           "<td style=\"border:2px solid #fff;\">success,</td>",
           "<td style=\"border:2px solid #fff;\"><div class=\"failure\"></div></td>",
           "<td style=\"border:2px solid #fff;\">failure,</td>",
           "<td style=\"border:2px solid #fff;\"><div class=\"error\"></div></td>",
           "<td style=\"border:2px solid #fff;\">error.</td></tr>",
           "</table><br>",
           ss, "")
  }
  
  
  ## 3: TAXONOMY
  ## ------------------
  ingroup.queried <- x@taxon@ingroup
  
  tax <- dbReadTaxonomy(x)
  ## the following code does not work with a parent-child table
  ## adapt or omit!
  # if (x@taxon@tip.rank == "gen"){
  #   id <- which(names(tax) == "gen")
  #   tax <- unique(tax[, 1:id])
  # } else {
  #   tax$spec <- gsub("_", " ", tax$spec)
  # }
  # taxonomy2html(tax)
  # ingroup.queried <- x@taxon@ingroup
  # if (spec.list){
  #   ingroup.queried <- sapply(ingroup.queried, head, n = 1)
  #   ingroup.received  <- tax$spec[grep("^ingroup", tax$tag)]
  #   extended.ingroup.received <- tax$spec[grep("extended ingroup", tax$tag)]
  #   ingroup.missing <- setdiff(ingroup.queried, ingroup.received)
  #   ingroup.false <- setdiff(ingroup.received, ingroup.queried)
  #   if ( length(ingroup.false) > 0 ){
  #     ingroup.false <- paste("- **WARNING**: *", ingroup.false, 
  #                            "* falsely as ingroup classified", sep = "")
  #   } else {
  #     ingroup.false <- ""
  #   }
  #   tax.coverage <- c(
  #     "Query of the NCBI Taxonomy Database gave these results:", "",
  #     paste("- **", length(ingroup.received), "** out of ", 
  #           length(ingroup.queried), " ingroup species (", 
  #           round(length(ingroup.received) / length(ingroup.queried) * 100, 2),
  #           "%) were found", sep = ""),
  #     paste("- **", length(ingroup.missing), "** ingroup species (", 
  #           round(length(ingroup.missing) / length(ingroup.queried) * 100, 2),
  #           "%) were missing", sep = ""),
  #     paste("- **", length(extended.ingroup.received), 
  #           "** congenerics of the ingroup species were included in the extended ingroup", 
  #           sep = ""),
  #     ingroup.false)
  #   
  # } else {
  #   tax.coverage <- c(
  #     "Query of the NCBI Taxonomy Database gave these results:", "",
  #     paste("- **", length(grep("ingroup", tax$tag)), 
  #           "** ingroup species", sep = ""),
  #     paste("- **", length(grep("outgroup", tax$tag)), 
  #           "** outgroup species", sep = "")
  #   )
  # }
  ## guide tree
  ## ----------
  # z <- c(z, "# Taxonomy",
  #        tax.coverage, "",
  #        "<a href='taxonomy.html'>view taxonomy table</a>",
  #        "",
  #        #          "```{r, echo=FALSE, message=FALSE}",
  #        #          "if ( inherits(x@taxon, 'taxonGuidetree') ){",
  #        #          "gt <- comprehensiveGuidetree(x, tip.rank = 'gen')",
  #        #          "} else {" ,
  #        #          "gt <- dbReadTaxonomy(x)",
  #        #          "gt <- tax2tree(gt, tip.rank = 'gen')",
  #        #          "}",
  #        #          "gt <- ladderize(gt)",
  #        #          "plot(gt, type = 'clado')",
  #        # "```", 
  #        "")
  
  ## stop here if no sequences have been downloaded so far
  ## -----------------------------------------------------
  if (!length(tabs)){
    write(z, file = file)
    return()  
  }
  
  ## RAW SEQUENCES
  ## -------------
  ll <- dbReadLocus(x)
  loci_genbank <- ll[, grep("gb_", names(ll))]
  loci_genbank[loci_genbank > 1] <- 1
  loci_genbank <- rowSums(loci_genbank)
  taxa_not_found <- names(loci_genbank)[loci_genbank == 0]
  n_taxa_not_found <- length(taxa_not_found)
  n_taxa_found <- nrow(ll) - n_taxa_not_found
  if (length(taxa_not_found)){
    taxaNotFound <- c(
      paste("##", Tip_rank ,"without any sequences"),
      paste0(n_taxa_not_found, " ", tip_rank, " (of ", nrow(ll), " available = ",
             round(100 * n_taxa_not_found/nrow(ll), 2), "%) have no sequences."))
    if (n_taxa_not_found <= 100){ 
      taxaNotFound <- c(taxaNotFound, paste0("- *", taxa_not_found, "*"))
    } 
  } else {
    taxaNotFound <- NULL
  }
  z <- c(z, "# DNA sequence retrieval",
         taxaNotFound, "",
         paste("## Number of", tip_rank, "per locus"),
         "```{r, echo=FALSE, message=FALSE}",
         "gb <- checkSpecLocus(x, 'gb')",
         "```", "")
  
  
  ## SELECTED SEQUENCES
  ## ------------------
  if (all(STATUS[, "F"] %in% c("pending", "failure", "error"))){
    write(z, file = file)
    return()  
  }
  loci_selected <- ll[, grep("sel_", names(ll))]
  loci_selected[loci_selected != "0"] <- 1
  loci_selected <- apply(loci_selected, c(1, 2), as.numeric) ## coerce to numeric!
  loci_selected <- rowSums(loci_selected)
  taxa_not_selected <- names(loci_selected)[loci_selected == 0]
  taxa_not_selected <- setdiff(taxa_not_selected, taxa_not_found)
  n_taxa_not_selected <- length(taxa_not_selected)
  if (length(taxa_not_selected)){
    taxaNotSelected <- c(
      paste("##", Tip_rank, "not selected for alignment"),
      paste0(n_taxa_not_selected, " ", tip_rank, " (of ", n_taxa_found, " retrieved = ",
             round(100* n_taxa_not_selected/n_taxa_found, 2), "%) have not been selected."))
    if (n_taxa_not_selected <= 100){
      taxaNotSelected <- c(taxaNotSelected, paste0("- *", taxa_not_selected, "*"))
    }
  } else {
    taxaNotSelected <- NULL
  }
  tab <- checkSpecLocus(x, "sel", plot = FALSE)
  tab <- cbind(rownames(tab$specPerMarker), tab$specPerMarker)
  colnames(tab) <- c("Locus", "Number of sequences", 
                     "Number of private sequences")
  tab <- htmlTable(tab) 
  
  z <- c(z, 
         "# Sequence selection and alignment",
         taxaNotSelected, " ",
         paste("## Number of selected", tip_rank, "per locus"),
         "```{r, echo=FALSE, message=FALSE}",
         "sel <- checkSpecLocus(x, 'sel')",
         "```", "", tab)
  
  if (all(STATUS[, "G"] %in% c("pending", "failure", "error"))){
    write(z, file = file)
    return()  
  }
  
  if (all(STATUS[, "H"] %in% c("pending", "failure", "error"))){
    write(z, file = file)
    return()  
  }
  
  ## SATURATION + BLOCKS
  ## ------------------
  msa <- dbSummaryMSA(x)
  msa[, -1] <- apply(msa[, -1, drop = FALSE], 2, 
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