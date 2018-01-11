## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-10-18)

#' @title Histograms of Identity/Coverage per Locus
#' 
#' @description Print histograms of the distribution of identity and 
#' coverage for each locus in the database to a PDF file.
#' @param megapteraProj An object of class \code{\link{megapteraProj}}.
#' @return None, as a side effect a PDF file is written to the working 
#' directory.
#' @seealso \code{\link{checkSpecLocus}},\code{\link{checkMissingSpec}}, 
#' \code{\link{checkSpecies}}.
#' @importFrom graphics hist par
#' @importFrom grDevices dev.off pdf
#' @export
 
checkIdentityCoverage <- function(megapteraProj){
  
  
  conn <- dbconnect(megapteraProj)
  acc <- DBI::dbListTables(conn)
  acc <- acc[grep("^acc_", acc)]
  
  pdf("identity-coverage.pdf")
  for ( i in acc ){
    SQL <- paste("SELECT identity, coverage",
                 "FROM", i)
    ic <- DBI::dbGetQuery(conn, SQL)
    par(mfrow = c(2, 1))
    hist(ic$identity, col = "yellow", xlab = "Identity",
         main = i)
    hist(ic$coverage, col = "skyblue", xlab = "Coverage", 
         main = NULL)
  }
  dev.off()
  system("open identity-coverage.pdf")
  
}

