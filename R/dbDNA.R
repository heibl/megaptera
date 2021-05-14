## This code is part of the megaptera package
## © C. Heibl 2017 (last update 2021-03-12)

#' @title Read and Write DNA Sequences
#' @description Read and write DNA Sequences from/to PostgreSQL Database.
#' @param x An object of class \code{\link{megapteraProj}}.
#' @param conn A connection object as produced by \code{\link[DBI]{dbConnect}}.
#' @param tab.name A vector of mode \code{"character"} giving the name of the
#'   database table. For \code{dbReadDNA} and when the locus is defined (see
#'   \code{\link{setLocus}}) this can be omitted and the data will be read from
#'   the species-level (or genus-level) database table.
#' @param taxon A vector of mode \code{"character"} used to choose a subset of
#'   available taxa. This can be either one or more taxon names or a regular
#'   expression.
#' @param regex Logical: if \code{TRUE}, the string given via \code{taxon} will
#'   be interpreted as a regular expression (see \code{\link{regex}}).
#' @param dna An object of class \code{\link[ape]{DNAbin}} to write to the
#'   database.
#' @param max.bp An integer, only sequences equal or shorter than \code{max.bp}
#'   will be returned.
#' @param max.evalue A positive real number giving the highest accepted Expectation
#'   value (E-value).
#' @param min.identity A real number between 0 and 1, only sequences with a
#'   fraction of at least \code{min.identity} nucleotides that are identical
#'   with the reference sequence will be returned.
#' @param min.coverage A real number between 0 and 1, only sequences with a
#'   fraction of at least \code{min.coverage} base pairs (compared to the
#'   reference sequence) will be returned.
#' @param status A vector of mode \code{"character"} to be written to the
#'   \emph{status} field of the PostgreSQL table.
#' @param ignore.excluded Logical: if \code{FALSE}, \code{dbReadDNA} will also
#'   return sequences that are marked in the \emph{status} field as 'excluded',
#'   'too long', or 'too distant'.
#' @param blocks A vector of mode \code{"character"} indicating how to handle
#'   alignment blocks: \code{"split"} causes blocks to be returned as elements
#'   of a list. \code{"concatenate"} means, blocks will be concatendated and
#'   returned as a single alignment. In order to get a (potentially unaligned)
#'   list of sequences, choose \code{"ignore"}.
#' @param reliability A vector of mode \code{"numeric"} between 0 and 1 giving
#'   the minimum reliability score for each alignment column.
#' @param subtree Logical.
#' @name dbDNA

NULL