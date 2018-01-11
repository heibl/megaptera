## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2018-01-11)

#' @title Class "dbPars"
#' @description S4 class for database parameters for a megaptera project, as
#'   defined with \code{\link{megapteraProj}}.
#' @slot host A vector of mode \code{"character"}, giving the name of the host,
#'   most often this will be \code{"localhost"}.
#' @slot port Numeric, giving the port number, most often \code{5432}.
#' @slot dbname A vector of mode \code{"character"}, giving the name of the
#'   database.
#' @slot user A vector of mode \code{"character"}, giving the name of the user.
#' @slot password A vector of mode \code{"character"}, giving the password.
#' @seealso \code{\link{dbPars}}, \code{\link{locus}}, and
#'   \code{\link{megapteraPars}} for defining of database parameters, loci, and
#'   the pipeline's parameters, respectively, and \code{\link{megapteraProj}}
#'   for the bundling of input data.

setClass("dbPars", 
         representation = list(
           host = "character", 
           port = "numeric", 
           dbname = "character", 
           user = "character", 
           password = "character")
)