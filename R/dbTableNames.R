## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2017-01-25)

#' @export
#' @import DBI

dbTableNames <- function(x, tabs = "all"){
  
  ## coerce full latin ranks to short form
  if (tabs == "species") tabs <- "spec"
  if (tabs == "genus") tabs <- "gen"

  conn <- dbconnect(x)
  tabnames <- paste("SELECT table_name",
                    "FROM information_schema.tables",
                    "WHERE table_schema='public'",
                    "AND table_type='BASE TABLE'")
  tabnames <- dbGetQuery(conn, tabnames)$table_name
  dbDisconnect(conn)
  
  ## select tables:
  tab.set <- c("all", "acc", x@taxon@tip.rank)
  tabs <- match.arg(tabs, tab.set)
  if (!"all" %in% tabs) {
    tabs  <- paste("^", tabs, "_", sep = "")
    tabnames <- tabnames[grep(tabs, tabnames)]
  }
  sort(tabnames)
}