## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-05-17)

#' @export
#' @import DBI

dbMaxGIPerSpec <- function(x, max.gi.per.spec){
  
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  logfile <- paste0("log/", gene, "-stepB.log")
  
  if (missing(max.gi.per.spec)) 
    max.gi.per.spec <- x@params@max.gi.per.spec
  
  conn <- dbconnect(x@db)
  taxon <- paste("SELECT gi, taxon, npos", 
                 "FROM", acc.tab,
                 "WHERE status != 'excluded (indet)'",
                 "AND status != 'excluded (too long)'")
  taxon <- dbGetQuery(conn, taxon)
  freqs <- table(taxon$taxon)
  if (any(id <- freqs > max.gi.per.spec)){
    id <- names(id)[id]
    slog("\n..", length(id), "species have >", 
         max.gi.per.spec, "sequenes:", 
         file = logfile)
    for (i in id){
      acc <- taxon[taxon$taxon == i, ]
      acc <- acc[order(acc$npos, decreasing = TRUE), ]
      acc.exclude <- tail(acc$gi, -max.gi.per.spec)
      acc.include <- setdiff(acc$gi, acc.exclude)
      slog("\n -", i, "(", length(acc.exclude), "sequences excluded )",
           file = logfile)
      acc <- paste("UPDATE", acc.tab, 
                   "SET status='raw'",
                   "WHERE", wrapSQL(acc.include, "gi", "=", "OR"))
      lapply(acc, dbSendQuery, conn = conn)
    }
  }
  dbDisconnect(conn)
}