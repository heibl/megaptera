## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2018-01-31)

#' @title RMarkdown Report
#' @description Creates a status report for a megaptera Project using RMarkdown.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param file A character vector giving a file name.
#' @param nmax An integer giving the maximun number of taxa for printed lists.
#' @return None, \code{megaptera2Rmarkdown} is called for its side effect.
#' @seealso \code{\link{checkSpecLocus}}, \code{\link{checkSpecies}}
#' @importFrom utils packageDescription
#' @export

megaptera2Rmarkdown <- function(x, file, nmax = 100){
  
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
  ingroup_speclist <- all(unlist(sapply(x@taxon@ingroup, is.Linnean)))
  ingroup <- head(x@taxon@ingroup, 3)
  if (ingroup_speclist) ingroup <- paste0("*", ingroup, "*")
  ingroup <- paste(ingroup, collapse = ", ")
  n <- length(x@taxon@ingroup)
  if (n > 3){
    ingroup <- paste0(ingroup, paste0(", ... [", n, " taxa in total]"))
  }
  outgroup_speclist <- all(unlist(sapply(x@taxon@outgroup, is.Linnean)))
  outgroup <- head(x@taxon@outgroup, 3)
  if (outgroup_speclist) outgroup <- paste0("*", outgroup, "*")
  outgroup <- paste(outgroup, collapse = ", ")
  n <- length(x@taxon@outgroup)
  if (n > 3){
    outgroup <- paste0(outgroup, paste0(", ... [", n, " taxa in total]"))
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
    names(STATUS) <- c("Locus", LETTERS[2:8])
    ss <- STATUS
    ss[ss == "pending"] <- 1
    ss[ss == "success"] <- 2
    ss[ss == "failure"] <- 3
    ss[ss == "error"] <- 4
    ss <- ss[order(ss[, 2], ss[, 3], ss[, 4], ss[, 5], ss[, 6], ss[, 7], ss[, 8]), ]
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
  tax <- dbReadTaxonomy(x)
  
  ## ingroup
  ## -------
  if (ingroup_speclist){
    ingroup_queried <- sapply(x@taxon@ingroup, head, n = 1)
    ingroup_received  <- intersect(ingroup_queried, tax$taxon)
    n_ingroup_received <- length(ingroup_received)
    ingroup_missing <- setdiff(ingroup_queried, ingroup_received)
    n_ingroup_missing <- length(ingroup_missing)
    ingroup_false <- setdiff(ingroup_received, ingroup_queried)
    
    ingroup_coverage <- c(
      "Query of the NCBI Taxonomy Database gave these results:", "",
      paste0("**", n_ingroup_received, "** ingroup species (of ",
             length(ingroup_queried), " = ",
             round(n_ingroup_received / length(ingroup_queried) * 100, 2),
             "%) could be retrieved"))
    if (n_ingroup_received <= nmax) ingroup_coverage <- c(ingroup_coverage, "", paste0("- *", ingroup_queried, "*"))
    if (n_ingroup_missing){
      ingroup_coverage <- c(ingroup_coverage, "",
                        paste0("**", n_ingroup_missing, "** ingroup species (",
                               round(n_ingroup_missing / length(ingroup_queried) * 100, 2),
                               "%) were missing"))
      if (n_ingroup_missing <= nmax) ingroup_coverage <- c(ingroup_coverage, "", paste0("- *", ingroup_missing, "*"))
      
    }
    if (length(ingroup_false)){
      ingroup_coverage <- c(ingroup_coverage, "",
                        paste0("- **WARNING**: *", ingroup_false, "* falsely as ingroup classified"))
    }
  } else {
    ranks <- data.frame(sg = c("species", "genus", "family", "order"),
                        pl = c("species", "genera", "families", "orders"),
                        stringsAsFactors = FALSE)
    ingroup <- taxdumpSubset(tax, x@taxon@ingroup, root = "mrca")
    ingroup_tab <- table(ingroup$rank)
    ig_ranks <- ranks[ranks$sg %in% names(ingroup_tab), ]
    ingroup_tab <- ingroup_tab[ig_ranks$sg]
    names(ingroup_tab)[ingroup_tab > 1] <- ig_ranks$pl[ingroup_tab > 1]
    ingroup_coverage <- c(paste0("NCBI Taxonomy was searched for **", x@taxon@tip.rank, 
                                  "** of **", x@taxon@ingroup, "**."), "",
                           paste0("The search retrieved ", 
                                  paste(paste(ingroup_tab, names(ingroup_tab)), collapse = " in ")), ".")
    if (ingroup_tab["species"] <= nmax) ingroup_coverage <- c(ingroup_coverage, "", 
                                                                paste0("- *", sort(ingroup$taxon[ingroup$rank == x@taxon@tip.rank]), "*"))
  }
  ingroup_coverage <- c("## Ingroup", ingroup_coverage)
  
  ## outgroup
  ## -------
  if (outgroup_speclist){
    outgroup_queried <- sapply(x@taxon@outgroup, head, n = 1)
    outgroup_received  <- intersect(outgroup_queried, tax$taxon)
    n_outgroup_received <- length(outgroup_received)
    outgroup_missing <- setdiff(outgroup_queried, outgroup_received)
    n_outgroup_missing <- length(outgroup_missing)
    outgroup_false <- setdiff(outgroup_received, outgroup_queried)
    
    outgroup_coverage <- c(
      "Query of the NCBI Taxonomy Database gave these results:", "",
      paste0("**", n_outgroup_received, "** outgroup species (of ",
             length(outgroup_queried), " = ",
             round(n_outgroup_received / length(outgroup_queried) * 100, 2),
             "%) could be retrieved"))
    if (n_outgroup_received <= nmax) outgroup_coverage <- c(outgroup_coverage, "", paste0("- *", outgroup_queried, "*"))
    if (n_outgroup_missing){
      outgroup_coverage <- c(outgroup_coverage, "",
                            paste0("**", n_outgroup_missing, "** outgroup species (",
                                   round(n_outgroup_missing / length(outgroup_queried) * 100, 2),
                                   "%) were missing"))
      if (n_outgroup_missing <= nmax) outgroup_coverage <- c(outgroup_coverage, "", paste0("- *", outgroup_missing, "*"))
      
    }
    if (length(outgroup_false)){
      outgroup_coverage <- c(outgroup_coverage, "",
                            paste0("- **WARNING**: *", outgroup_false, "* falsely as outgroup classified"))
    }
  } else {
    ranks <- data.frame(sg = c("species", "genus", "family", "order"),
                        pl = c("species", "genera", "families", "orders"),
                        stringsAsFactors = FALSE)
    outgroup <- taxdumpSubset(tax, x@taxon@outgroup, root = "mrca")
    outgroup_tab <- table(outgroup$rank)
    og_ranks <- ranks[ranks$sg %in% names(outgroup_tab), ]
    outgroup_tab <- outgroup_tab[og_ranks$sg]
    names(outgroup_tab)[outgroup_tab > 1] <- og_ranks$pl[outgroup_tab > 1]
    outgroup_coverage <- c(paste0("NCBI Taxonomy was searched for **", x@taxon@tip.rank, 
                                 "** of **", x@taxon@outgroup, "**."), "",
                          paste0("The search retrieved ", 
                                paste(paste(outgroup_tab, names(outgroup_tab)), collapse = " in "), "."))
    if (outgroup_tab["species"] <= nmax) outgroup_coverage <- c(outgroup_coverage, "", 
                                                                paste0("- *", sort(outgroup$taxon[outgroup$rank == x@taxon@tip.rank]), "*"))
  }
  outgroup_coverage <- c("## Outgroup", outgroup_coverage)
  
  ## guide tree
  ## ----------
  if (inherits(x@taxon, 'taxonGuidetree')){
    guide_tree <- c("## User-defined guide tree",
                    "```{r, echo=FALSE, message=FALSE}",
                    "plot(x@taxon@guide.tree, type = 'phylo', no.margin = TRUE)",
                    "```", "")
  } else {
    guide_tree <- ""
  }
  z <- c(z, "# Taxonomy", ingroup_coverage, "", outgroup_coverage, "", guide_tree)
  
  
  
  ## stop here if no sequences have been downloaded so far
  ## -----------------------------------------------------
  if (!length(tabs)){
    write(z, file = file)
    return()  
  }
  
  ## RAW SEQUENCES
  ## -------------
  ll <- dbReadLocus(x)
  loci_genbank <- ll[grep("gb", names(ll))]
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
    if (n_taxa_not_found <= nmax){ 
      taxaNotFound <- c(taxaNotFound, "", paste0("- *", taxa_not_found, "*"))
    } 
  } else {
    taxaNotFound <- NULL
  }
  z <- c(z, "# DNA sequence retrieval",
         taxaNotFound, "",
         paste("## Number of", tip_rank, "per locus"),
         "```{r, echo=FALSE, message=FALSE}",
         "gb <- checkSpecLocus(x, 'retrieved')",
         "```", "")
  
  
  ## SELECTED SEQUENCES
  ## ------------------
  if (all(STATUS[, "F"] %in% c("pending", "failure", "error"))){
    write(z, file = file)
    return()  
  }
  loci_selected <- ll[, grep("sel_", names(ll)), drop = FALSE]
  loci_selected[loci_selected != "0"] <- 1
  loci_selected <- apply(loci_selected, c(1, 2), as.numeric) ## coerce to numeric!
  loci_selected <- rowSums(loci_selected)
  taxa_not_selected <- names(loci_selected)[loci_selected == 0]
  taxa_not_selected <- setdiff(taxa_not_selected, taxa_not_found)
  n_taxa_not_selected <- length(taxa_not_selected)
  if (n_taxa_not_selected){
    taxaNotSelected <- c(
      paste("##", Tip_rank, "not selected for alignment"),
      paste0(n_taxa_not_selected, " ", tip_rank, " (of ", nrow(ll), " available = ", 
             round(100* n_taxa_not_selected/nrow(ll), 2),
             "% and of ", n_taxa_found, " retrieved = ",
             round(100* n_taxa_not_selected/n_taxa_found, 2), "%) have not been selected."))
    if (n_taxa_not_selected <= nmax){
      taxaNotSelected <- c(taxaNotSelected, "", paste0("- *", taxa_not_selected, "*"))
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