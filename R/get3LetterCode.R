## This code is part of the megaptera package
## © C. Heibl 2019 (last update 2019-11-13)

#' @title GenBank Division Code
#' @description Get the GenBank Division Code based on the current taxonomy.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @return A vector of mode \code{"character"} containing three-letter
#'   abreviations of the GenBank Division Code.
#' @references \url{https://www.ncbi.nlm.nih.gov/genbank/htgs/divisions/}
#' @export

get3LetterCode <- function(x){
  
  tax <- dbReadTaxonomy(x, root = "tol")
  out <- NULL
  
  ## Mammals
  ## -------
  if ("Mammalia" %in% tax$taxon){
    o <- taxdumpChildren(tax, "Mammalia", tip.rank = "order")
    o <- o$taxon[o$rank == "order"]
    if ("Primates" %in% o){
      out <- c(out, "PRI")
    }
    if ("Rodentia" %in% o){
      out <- c(out, "ROD")
    }
    if (any(!o %in% c("Primates", "Rodentia")))
      out <- c(out, "MAM")
  }
  
  ## Vertebrates other than Mammals
  ## ------------------------------
  if ("Vertebrata" %in% tax$taxon){
    m <- taxdumpMRCA(tax) ## actual root
    if (!"Mammalia" %in% taxdumpLineage(tax, m)$taxon){
      out <- c(out, "VRT")
    }
  }
  
  ## Invertebrates
  ## ------------------------------
  if ("Metazoa" %in% tax$taxon){
    m <- taxdumpMRCA(tax) ## actual root
    if (!"Vertebrata" %in% taxdumpLineage(tax, m)$taxon){
      out <- c(out, "INV")
    }
  }
  
  ## Plants, algi and fungi
  ## ----------------------
  tt <- c("Fungi", "Rhodophyta", "Stramenopiles", "Viridiplantae")
  if (any(tt %in% tax$taxon)){
      out <- c(out, "PLN")
  }
  
  # code <- rbind(
  #   c("PRI", "primate sequences"),
  #   c("ROD", "rodent sequences"),  
  #   c("MAM", "other mammalian sequences"),
  #   c("VRT", "other vertebrate sequences"),
  #   c("INV", "invertebrate sequences"),
  #   c("PLN", "plant, fungal, and algal sequences"),
  #   c("BCT", "bacterial sequences"),
  #   c("VRL", "viral sequences"),
  #   c("PHG", "bacteriophage sequences"),
  #   c("SYN", "synthetic sequences"),
  #   c("UNA", "unannotated sequences"),
  #   c("EST", "EST sequences (Expressed Sequence Tags)"),
  #   c("PAT", "patent sequences"),
  #   c("STS", "STS sequences (Sequence Tagged Sites)"),
  #   c("GSS", "GSS sequences (Genome Survey Sequences)"),
  #   c("HTG", "HTGS sequences (High Throughput Genomic sequences)"),
  #   c("HTC", "HTC sequences (High Throughput cDNA sequences)"),
  #   c("ENV", "Environmental sampling sequences"),
  #   c("CON", "Constructed sequences"),
  #   c("TSA", "Transcriptome Shotgun Assembly sequences"))
  out
}