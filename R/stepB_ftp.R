## This code is part of the megaptera package
## © C. Heibl 2019 (last update 2019-09-17)

#' @import DBI
#' @importFrom restez db_download restez_path_set
#' @export

stepB_ftp <- function(x, update.seqs = "no"){
  
  
 cat("Switch works!") 
  
  restez_path_set("data")
  
  db_download #(restez::db-setup-tools)
}