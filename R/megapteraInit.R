## This code is part of the megaptera package
## © C. Heibl 2019 (last update 2021-03-12)

#' @title Initialize MEGAPTERA
#' @description Creates a consist system of folders where data and metadata are
#'   stored.
#' @param path A vector of mode \code{"character"} that gives the directory
#'   where megaptera should store data and results.
#' @details \code{"megapteraInit"} creates a file system that is designed to
#'   keep memory requirements and execution times at the necessary minimum
#'   during subsequent runs of the MEGAPTERA pipeline. Individual projects can
#'   be created with \code{\link{megapteraProj}}.
#' @return None, \code{megapteraInit} is called for its side effects.
#' @seealso \code{\link{megapteraProj}}, \code{\link{ncbiTaxonomy}}
#' @import DBI
#' @importFrom restez db_download restez_path_set
#' @export

megapteraInit <- function(path){
  
  megaptera_path <- file.path(path, "megaptera_data")
  taxonomy_path <- file.path(megaptera_path, "gb_taxonomy")
  sequence_path <- file.path(megaptera_path, "gb_sequence")
  project_path <- file.path(megaptera_path, "project")
  
  ## root directory 'megaptera_data'
  ## -------------------------------
  if (!dir.exists(megaptera_path)){
    dir.create(megaptera_path)
    info <- c(as.character(Sys.time()),
              packageDescription("megaptera")$Version)
    info <- paste(info, collapse = "\n")
    write(info, file.path(megaptera_path, "INFO"))
  } else {
    info <- scan(file.path(megaptera_path, "INFO"), 
                 what = "c", sep = "\n", quiet = TRUE)
      message("megaptera file system was already initialized on ",
              info[1], " by package version ", info[2])
  }
  
  ## Make sure that subdirectories exist
  ## -----------------------------------
  if (!dir.exists(taxonomy_path)){
    dir.create(taxonomy_path)
  }
  if (!dir.exists(sequence_path)){
    dir.create(sequence_path)
  }
  if (!dir.exists(project_path)){
    dir.create(project_path)
  }
}