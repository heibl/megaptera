## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-05-11)

#' @title Limit Numbers of Sequences per Species
#' @description Some species (e.g. model organism) have thousands or more
#'   sequences of a single locus on GenBank. Much of this genetic information is
#'   redundant in a phylogenetic context, so \strong{megaptera} internally
#'   limits the number of sequences per species to 10. This number can be
#'   changed by the user via with \code{megapteraPars} or \code{dbMaxGIPerSpec},
#'   the latter function allowing more fine-tuning (see Details).
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @param max.gi.per.spec An integer giving the maximum number of sequences per
#'   species to be used by the pipeline.
#' @param prefer A character string indicating under what criterion the
#'   sequences will be choosen; available are \code{"longest"} (default),
#'   \code{"shortest"}, \code{"most frequent length"} and \code{"random"}.
#' @param taxon A character string giving one or more taxon names for which the
#'   number of sequences will be limited. If left missing, all taxon names found
#'   in the database are handled.
#' @details After calling \code{dbMaxGIPerSpec} the pipeline must be rerun from
#'   \code{\link{stepC}} onwards.
#'
#'   \emph{To do: describe parameters!}
#' @seealso \code{\link{megapteraPars}}
#' @import DBI
#' @importFrom utils tail
#' @export


dbMaxGIPerSpec <- function(megProj, max.gi.per.spec, prefer = "longest", taxon){
  
  ## PARAMETERS
  ## -----------
  gene <- megProj@locus@sql
  if (gene == "undefined") stop("undefined locus not allowed")
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  logfile <- paste0("log/", gene, "-stepB.log")
  prefer <- match.arg(prefer, c("longest", "shortest", "most frequent length", "random"))
  
  if (missing(max.gi.per.spec)) 
    max.gi.per.spec <- megProj@params@max.gi.per.spec
  
  conn <- dbconnect(megProj@db)
  
  ## Delete entries from previous runs
  ## ---------------------------------
  SQL <- paste("UPDATE", acc.tab,
               "SET status='raw'",
               "WHERE status='excluded (max.gi)'")
  if (!missing(taxon)) {
    SQL <- paste(SQL, "AND", wrapSQL(gsub("_", " ", taxon), "taxon", "="))
  }
  dbSendQuery(conn, SQL)
  
  
  ## Exclude accession (if necessary)
  ## --------------------------------
  tax <- paste("SELECT gi, taxon, npos", 
                 "FROM", acc.tab,
                 "WHERE status != 'excluded (indet)'",
                 "AND status != 'excluded (too long)'")
  
  ## handle only a subset of taxa
  ## ----------------------------
  if (!missing(taxon)) {
    tax <- paste(tax, "AND", wrapSQL(gsub("_", " ", taxon), "taxon", "="))
  }
  
  tax <- dbGetQuery(conn, tax)
  freqs <- table(tax$taxon)
  if (any(id <- freqs > max.gi.per.spec)){
    id <- names(id)[id]
    if (missing(taxon)){
      slog("\n..", length(id), "species have >", 
           max.gi.per.spec, "sequenes:", 
           file = logfile)
    }
    for (i in id){
      acc <- tax[tax$taxon == i, ]
      if (prefer == "longest"){
        acc <- acc[order(acc$npos, decreasing = TRUE), ]
        acc_exclude <- tail(acc$gi, -max.gi.per.spec)
      }
      if (prefer == "shortest"){
        acc <- acc[order(acc$npos, decreasing = FALSE), ]
        acc_exclude <- tail(acc$gi, -max.gi.per.spec)
      }
      if (prefer == "most frequent length"){
        freqs <- sort(table(acc$npos), decreasing = TRUE)
        freqs[] <- 1:length(freqs)
        acc$order <- freqs[match(acc$npos, names(freqs))]
        acc <- acc[order(acc$order), ]
        acc_exclude <- tail(acc$gi, -max.gi.per.spec)
      }
      if (prefer == "random"){
        acc_exclude <- sample(acc$gi, max.gi.per.spec)
      }
      acc_include <- setdiff(acc$gi, acc_exclude)
      
      ## prepare and execute SQL statements
      ## ----------------------------------
      slog("\n -", i, "(", length(acc_exclude), "sequences excluded )",
           file = logfile)
      acc_include <- paste("UPDATE", acc.tab,
                           "SET status='raw'",
                           "WHERE", wrapSQL(acc_include, "gi", "=", "OR"))
      lapply(acc_include, dbSendQuery, conn = conn)
      acc_exclude <- paste("UPDATE", acc.tab,
                           "SET status='excluded (max.gi)', identity=NULL, coverage=NULL",
                           "WHERE", wrapSQL(acc_exclude, "gi", "=", "OR"))
      lapply(acc_exclude, dbSendQuery, conn = conn)
      
    }
  }
  invisible(dbDisconnect(conn))
}