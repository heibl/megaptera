## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-03-07)

#' @export

dbProgress <- function(megProj, step, status, locus.wise = TRUE){
  
  conn <- dbconnect(megProj)
  gene <- megProj@locus@sql
  SQL <- "SELECT * FROM progress"
  if (locus.wise & x@locus@sql != "undefined") {
    SQL <- paste(SQL, "WHERE", wrapSQL(gene, "locus", "="))
  }
  y <- dbGetQuery(conn, SQL)
  
  if (!missing(status)){
    if (nrow(y)){
      SQL <- paste("UPDATE progress",
                   "SET", wrapSQL(status, step, "="), 
                   "WHERE", wrapSQL(gene, "locus", "="))
    } else {
      SQL <- wrapSQL(c(gene, "error", rep("pending", 7)), NULL, boolean = ",")
      SQL <- paste("INSERT INTO progress", 
                   "(locus, step_b, step_c, step_d, step_e, step_f, step_g, step_h, step_i)",
                   "VALUES (", SQL, ")")
  }
    dbSendQuery(conn, SQL)
  }
  dbDisconnect(conn)
  y
}