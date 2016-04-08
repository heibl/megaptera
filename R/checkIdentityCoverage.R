## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-05-07)

checkIdentityCoverage <- function(megapteraProj){
  
  
  conn <- dbconnect(megapteraProj)
  acc <- dbListTables(conn)
  acc <- acc[grep("^acc_", acc)]
  
  pdf("identity-coverage.pdf")
  for ( i in acc ){
    SQL <- paste("SELECT identity, coverage",
                 "FROM", i)
    ic <- dbGetQuery(conn, SQL)
    par(mfrow = c(2, 1))
    hist(ic$identity, col = "yellow", xlab = "Identity",
         main = i)
    hist(ic$coverage, col = "skyblue", xlab = "Coverage", 
         main = NULL)
  }
  dev.off()
  system("open identity-coverage.pdf")
  
}

