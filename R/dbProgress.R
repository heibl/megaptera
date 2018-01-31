## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-01-31)

#' @importFrom DBI dbDisconnect dbGetQuery dbSendQuery
#' @export

dbProgress <- function(megProj, step, status = "success", locus.wise = TRUE){
  
  conn <- dbconnect(megProj)
  gene <- megProj@locus@sql
  SQL <- "SELECT * FROM progress"
  if (locus.wise & gene != "undefined") {
    SQL <- paste(SQL, "WHERE", wrapSQL(gene, "locus", "="))
  }
  y <- dbGetQuery(conn, SQL)
  
  ## Update relation 'progress'
  ## --------------------------
  if (!missing(step)){
    if (gene == "undefined") stop("undefined locus not allowed")
    if (nrow(y)){
      ## Subsequent entries into relation 'progress'
      ## -------------------------------------------
      SQL <- paste("UPDATE progress",
                   "SET", wrapSQL(status, step, "="), 
                   "WHERE", wrapSQL(gene, "locus", "="))
    } else {
      ## First entry of locus into relation 'progress' (by definition stepB)
      ## -------------------------------------------------------------------
      SQL <- wrapSQL(c(gene, status, rep("pending", 6)), NULL, boolean = ",")
      SQL <- paste("INSERT INTO progress", 
                   "(locus, step_b, step_c, step_d, step_e, step_f, step_g, step_h)",
                   "VALUES (", SQL, ")")
    }
    dbSendQuery(conn, SQL)
  }
  
  dbDisconnect(conn)
  y
}