## This code is part of the megaptera package
## © C. Heibl 2015 (last update 2015-12-23)

dbTableNames <- function(x, tabs = "all"){

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
  if ( !"all" %in% tabs ) {
    tabs  <- paste("^", tabs, "_", sep = "")
    tabnames <- tabnames[grep(tabs, tabnames)]
  }
  sort(tabnames)
}