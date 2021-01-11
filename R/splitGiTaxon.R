## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2020-02-23)

#' @export

splitGiTaxon <- function(x, enforce.binomial = FALSE, sep = " "){
  
  # old version until 2020-02-23
  # fun <- function(x){
  #   x <- unlist(strsplit(x, "_"))
  #   c(gi = tail(x, 1), 
  #     taxon = paste(head(x, -1), collapse = sep))
  # }
  # x <- lapply(x, fun)
  # x <- do.call(rbind, x)
  # x <- as.data.frame(x, stringsAsFactors = FALSE)
  # 
  # if (enforce.binomial){
  #   x[, 2] <- strip.infraspec(x[, 2])
  # }
  
  reg_exp <- "(^[[:upper:]][[:lower:][:space:]_-]+)([ _])([[:upper:][:digit:]_ ]+$)"
  list(taxon = gsub(reg_exp, "\\1", x), gi = gsub(reg_exp, "\\3", x))
}
