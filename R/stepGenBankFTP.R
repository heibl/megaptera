## This code is part of the megaptera package
## © C. Heibl 2019 (last update 2019-11-13)

#' @title Extract Records from GenBank Flatfiles
#' @description Extract sequences of study taxa from GenBank flatfile.
#' @param x An object of class \code{\link{megaperaProj}}
#' @param update A vector of mode \code{"character"} giving the name of a taxon.
#' @import DBI
#' @importFrom ips blastn
#' @importFrom parallel clusterEvalQ clusterExport makeCluster parLapply
#'   parSapply stopCluster
#' @importFrom restez db_download restez_path_set
#' @export

stepGenBankFTP <- function(x, update){
  
  # start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  
  ## List of taxa to extract
  ## -----------------------
  taxa <- c(x@taxon@ingroup, x@taxon@outgroup)
  
  code <- get3LetterCode(x)
  # if (length(code) > 1) stop("debug me!")
  message("Division code: ", code)
  
  ## Hier weitermachen: Files herunterladen/ Welche Files müssen noch heruntergeladen werden
  ## code in : identify_downloadable_files(), nicht exportiert

  release <- readLines(url('ftp://ftp.ncbi.nlm.nih.gov/genbank/gbrel.txt'))
  
  file_size <- grep("File Size[[:space:]]+File Name", release)
  id <- grep("gbvrt[[:digit:]]+[.]seq", release)
  id <- id[id > file_size]
  seq_files <- gsub("^ +", "", release[id])
  seq_files <- do.call(rbind, strsplit(seq_files, "[[:space:]]+"))
  seq_files <- as.data.frame(seq_files, stringsAsFactors = FALSE)
  names(seq_files) <- c("size", "name")
  
  path <- file.path("/Volumes/HD710 PRO", "megaptera_data/gb_sequence")
  
  for (i in seq_files$name[2:nrow(seq_files)]){
    download.file(file.path("ftp://ftp.ncbi.nlm.nih.gov/genbank", paste0(i, ".gz")),
                  file.path(path, paste0(i, ".gz")))
  }
  
  
  
  
  lapply(list.files(path = "/Users/heibl/Documents/r/pkgs/restez-master/R", full.names = TRUE), source)
  restez_path_set(filepath = file.path(x@params@data.path, "megaptera_data/gb_sequence"))
  db_download(db = 'nucleotide', preselection = code[1])
  # Interactively download GenBank data 
  
  
  taxa <- dbReadTaxonomy(x)
  taxa <- taxa$taxon[taxa$rank == x@taxon@tip.rank & taxa$status == "scientific name"]

  ## 2. Read flat file into R based on restez::flatfile_read
  ##    either sequential or parallel
  ## -------------------------------------------------------
  # slog()
  
  rec <- paste0("gb", tolower(code), "[[:digit:]]+.seq.gz")
  rec <- paste(rec, collapse = "|")
  rec <- list.files(path = "../../gb_sequence/", ## "../../gb_sequence/restez/downloads/"
                    pattern = rec, full.names = TRUE)
  message("Found ", length(rec), " GenBank flatfiles")
  # file <- rec <- rec[grep("gbpln97[.]", rec)]
  # id <- sapply(rec, findInFlatfile, acc = "M32501")
  
  ## 3. In update modus identify flatfiles that contain the taxon
  ##    that should be updated
  ## -------------------------
  if (!missing(update)){
    taxa <- taxa[grep(update, taxa)]
    cl <- makeCluster(x@params@cpus)
    clusterEvalQ(cl, library(megaptera))
    clusterExport(cl, c("taxa", "rec"), environment())
    id <- parSapply(cl, X = rec, FUN = findInFlatfile, taxa = unlist(taxa)) # 12:13
    stopCluster(cl)
    
    rec <- names(id)[id]
  }
  
  ## 4. Extract sequences from flatfiles and write them to database
  ## --------------------------------------------------------------
  if (length(rec) < x@params@cpus | !x@params@parallel){
    lapply(rec, gbflat2db, taxa = taxa, megProj = x)
  } else {
    cl <- makeCluster(x@params@cpus)
    clusterEvalQ(cl, library(megaptera))
    clusterExport(cl, c("x", "rec", "taxa"), environment())
    parLapply(cl, X = rec, fun = gbflat2db, taxa = taxa, megProj = x)
    stopCluster(cl)
  }
}
