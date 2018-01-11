## This code is part of the megaptera package.
## It is based on code of Scott Scamberlein and
## streamlined towards the needs of the megaptera 
## package
## © C. Heibl 2017 (last update 2017-10-23)

#' @title Barcode of Life Database Interface
#' @description Get DNA sequences from BOLD SYSTEMS
#' @param taxon A character string giving a taxon name.
#' @param marker A character string giving a marker name.
#' @export

BOLD2megaptera <- function (taxon, marker) {
  
  args <- list(taxon = taxon, marker = marker)
  cli <- crul::HttpClient$new("http://v4.boldsystems.org/index.php/API_Public/sequence")
  out <- cli$get(query = args)
  out$raise_for_status()
  if (grepl("html", out$response_headers$`content-type`)) {
    stop(out$parse("UTF-8"))
  }
  out <- out$parse("UTF-8")
  out <- strsplit(out, ">")[[1]][-1]
  
  
  out <- gsub("([\n])([[:upper:]-])", "|", out)
  out <- gsub("[\r\n]", "", out)
  out <- strsplit(out, "|", fixed = TRUE)
  fiver <- function(z) z[c(1, 2, 3, NA, 4)]
  id <- sapply(out, length) == 4
  out[id] <- lapply(out[id], fiver)
  out <- do.call(rbind, out)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  names(out) <- c("id", "taxon", "marker", "ncbi", "sequence")
  out
}

  
  
  