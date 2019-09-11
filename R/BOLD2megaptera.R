## This code is part of the megaptera package.
## It is based on code of Scott Scamberlein and
## streamlined towards the needs of the megaptera 
## package
## © C. Heibl 2017 (last update 2019-04-16)

#' @title Barcode of Life Database Interface
#' @description Get DNA sequences from BOLD SYSTEMS
#' @param taxon A character string giving a taxon name (only one allowed).
#' @param marker A character string giving a marker name.
#' @param out.format A character string giving the desired output format;
#'   possible are \code{"df"} (data frame) or \code{"\link{DNAbin}"}
#' @export

BOLD2megaptera <- function (taxon, marker, out.format = "df") {
  
  args <- list(taxon = taxon, marker = marker)
  cli <- crul::HttpClient$new("http://v4.boldsystems.org/index.php/API_Public/sequence")
  out <- cli$get(query = args)
  out$raise_for_status()
  if (grepl("html", out$response_headers$`content-type`)) {
    stop(out$parse("UTF-8"))
  }
  out <- out$parse("UTF-8")
  out <- strsplit(out, ">")[[1]][-1]
  
  ## Construct data frame
  ## --------------------
  out <- gsub("([\n])([[:upper:]-])", "|", out)
  out <- gsub("[\r\n]", "", out)
  out <- strsplit(out, "|", fixed = TRUE)
  fiver <- function(z) z[c(1, 2, 3, NA, 4)]
  id <- sapply(out, length) == 4
  out[id] <- lapply(out[id], fiver)
  out <- do.call(rbind, out)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  names(out) <- c("id", "taxon", "marker", "ncbi", "sequence")
  
  ## Convert to class DNAbin if desired
  ## ----------------------------------
  if (out.format == "DNAbin"){
    bb <- paste(out[, "taxon"], out[, "id"])
    bb <- gsub(" ", "_", bb)
    out <- as.list(out[, "sequence"])
    out <- lapply(out, strsplit, split = "")
    out <- lapply(out, unlist)
    out <- lapply(out, tolower)
    names(out) <- bb
    out <- as.DNAbin(out)
  }
  
  out
}

  
  
  