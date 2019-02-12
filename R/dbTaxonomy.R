## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2017-10-23)

#' @title Read and Write Taxonomy Table
#' @description Create, extend and read the \code{taxonomy} table of a
#'   \code{megaptera} project.
#' @param megProj An object of class \code{\link{megapteraProj}}.
#' @param taxonomy A data frame containing a taxonomic classification of the
#'   study species. It will \bold{replace any existing} taxonomy in the database.
#' @param tag A character string used to tag the rows of \code{taxonomy} in the
#'   postgreSQL table; this argument is mostly for internal use.
#' @param tip.rank A character string giving the rank of the tips (e.g.
#'   \code{species}, \code{genus}, ...).
#' @param subset A subset of species names (Latin binomials) to which the
#'   taxonomy should be limited. Can be a DNA alignment of class \code{DNAbin},
#'   a phylogenetic tree of class \code{phylo},\code{species_sequence} (see
#'   Details), or a character vector.
#' @param root A character string. By default (\code{"tol"}),  all higher ranks
#'   from the most recent common ancestor (MRCA) of the clade of interest down
#'   to the root are included in the taxonomic classification, whereas
#'   \code{"mrca"} means that only the MRCA will be included. Note that this
#'   option only makes sense in combination with \code{subset}.
#' @param logfile A character string giving the names of a potential logfile; no
#'   logfile is written if this argument is left empty.
#' @return For \code{dbReadTable} an object of class \code{"data.frame"} holding
#'   the taxonomic classification of target species. \code{dbUpdateTaxonomy} is
#'   called only for its side effect creating and extending the \code{taxonomy}
#'   database table.
#'
#' @details If step \code{\link{stepF}} has been run, \code{dbReadTaxonomy} can
#'   be called with \code{subset = "species_sequence"}. In that case all species
#'   names with sequence data for the current locus will be used for subsetting.
#' @seealso \code{\link{ncbiTaxonomy}} for retrieval of taxonomies from the
#'   taxonomy database at NCBI; \code{\link{findRoot}} returns the lowest common
#'   rank for taxonomies. \code{\link{dbReadDNA}} and \code{\link{dbWriteDNA}}
#'   to read and write DNA sequences.
#' @name dbTaxonomy


NULL
