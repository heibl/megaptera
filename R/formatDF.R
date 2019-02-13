## This code is part of the megaptera package
## © C. Heibl 2018 (last update 2018-12-19)

#' @export

## Prepare dataframe for cat()
## ---------------------------
formatDF <- function(df){
  
  df <- rbind(names(df), df)
  df <- apply(df, 2, format, justify = "r")
  df <- apply(df, 1, paste, collapse = " | ")
  df <- c(df[1], paste(rep("-", nchar(df[1])), collapse = ""), df[2:length(df)])
  paste0("\n", df)
}