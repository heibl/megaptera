## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-01-30)

#' @title Read Sequences from Database
#' @description Reads selected and assembled sequences (see \code{\link{stepF}}
#'   and subsequent steps) from the database, specifically from the table
#'   'species_sequences'.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param taxon A vector of mode \code{"character"} used to choose a subset of
#'   available taxa. This can be either one or more taxon names or a regular
#'   expression.
#' @param regex Logical: if \code{TRUE}, the string given via \code{taxon} will
#'   be interpreted as a regular expression (see \code{\link{regex}}).
#' @param ignore.excluded \emph{Currently unused}.
#' @param blocks \emph{Currently unused}.
#' @return An object of class \code{\link{DNAbin}}.
#' @seealso \code{\link{dbReadDNA}} and \code{\link{dbWriteDNA}}
#' @importFrom ape as.DNAbin cbind.DNAbin
#' @importFrom DBI dbDisconnect dbGetQuery
#' @export

dbReadMSA <- function(x, taxon, regex = TRUE, 
                      ignore.excluded = TRUE, 
                      blocks = "ignore"){
  
  ## Check and prepare input data
  ## ----------------------------
  if (!inherits(x, "megapteraProj")) stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("locus undefined; use setLocus() to define a locus")
  gene <- x@locus@sql
  tip.rank <-x@taxon@tip.rank
  msa.tab <- paste(tip.rank, "sequence", sep = "_")
  
  
  if (missing(taxon)) taxon <- ".+"
  otaxon <- taxon
  blocks <- match.arg(blocks, c("ignore", "split", "concatenate"))
  
  ## Escape metacharacters in taxon names
  ## -------------------------------------
  if (!regex){
    taxon <- gsub("_", " ", taxon)
    taxon <- gsub("[.]$", "[.]", taxon)
    taxon <- gsub("([(]|[+]|[)])", "[\\1]", taxon)
    taxon <- gsub("'", ".", taxon) # e.g.Acorus_sp._'Mt._Emei'
    taxon <- paste0("^", taxon, "$")
    taxon <- paste(taxon, collapse = "|")
  }
  
  ## Open database connection
  conn <- dbconnect(x)
  
  ## retrieve sequences from species.table or genus.table
  ## ----------------------------------------------------
  SQL <- paste("SELECT regexp_replace(taxon, ' ', '_') AS taxon, sequence, reliability",
               "FROM", msa.tab,
               "WHERE", wrapSQL(gene, "locus", "="))
  seqs <- dbGetQuery(conn, SQL)
  dbDisconnect(conn)
  if (!nrow(seqs)) {
    warning("no sequences for\n- ", paste(otaxon, collapse = "\n- "))
    return(NULL)
  }
  seqs <- pg2DNAbin(seqs)
  
  # ## split into blocks (if necessary)
  # ## --------------------------------
  # if (length(grep("block 2", block)) & blocks %in% c("split", "concatenate")){
  #   
  #   block <- split(block$taxon, f = block$status)
  #   seqs <- lapply(block, 
  #                  function(a, tips) deleteEmptyCells(as.matrix(a[tips]), 
  #                                                     quiet = TRUE),
  #                  a = seqs)
  #   if (blocks == "concatenate"){
  #     seqs <- do.call(cbind.DNAbin, c(seqs, fill.with.gaps = TRUE))
  #   }
  # } else {
  #   
  #   ## if single block: convert to matrix if possible
  #   ## ----------------------------------------------
  #   if (length(unique(sapply(seqs, length))) == 1){
  #     seqs <- as.matrix(seqs)
  #     seqs <- deleteEmptyCells(seqs, quiet = TRUE)
  #   }
  # }
  return(seqs)  
}
