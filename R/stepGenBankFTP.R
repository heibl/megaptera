## This code is part of the megaptera package
## © C. Heibl 2019 (last update 2021-03-17)

#' @title Extract Records from GenBank Flatfiles
#' @description Extract sequences of study taxa from GenBank flatfile.
#' @param x An object of class \code{\link{megapteraProj}}
#' @param update A vector of mode \code{"character"} giving the name of a taxon.
#' @import DBI
#' @importFrom crayon %+% bold cyan green magenta silver
#' @importFrom ips blastn
#' @importFrom lubridate mdy
#' @importFrom parallel clusterEvalQ clusterExport makeCluster parLapply
#'   parSapply stopCluster
#' @importFrom restez db_download restez_path_set
#' @export

stepGenBankFTP <- function(x, update){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  
  ## iniate logfile
  ## --------------
  # logfile <- paste0("log/", gene, "-stepF.log")
  # if ( file.exists(logfile) ) unlink(logfile)
  slog(silver(bold(paste("megaptera", packageDescription("megaptera")$Version)) %+% "\n"
              %+% as.character(Sys.time()) %+%  "\n"
              %+% bold("STEP GenBankFTP") %+% ": download sequences flat files from GenBabnk FTP server\n"))
  
  ## List of taxa to extract
  ## -----------------------
  taxa <- c(x@taxon@ingroup, x@taxon@outgroup)
  
  code <- get3LetterCode(x)
  # if (length(code) > 1) stop("debug me!")
  cat(silver("Division code: " %+% magenta$bold(code) %+% "\n"))
  
  ## Hier weitermachen: Files herunterladen/ Welche Files müssen noch heruntergeladen werden
  ## code in : identify_downloadable_files(), nicht exportiert
  options(timeout = max(300, getOption("timeout"))) ## default of 60 s is too short
  conn <- url('ftp://ftp.ncbi.nlm.nih.gov/genbank/gbrel.txt')
  release <- readLines(conn)
  close(conn)
  
  ## Number and date of latest release
  ## ------------------------------------
  rel_nr <- grep("NCBI-GenBank Flat File Release", head(release))
  rel_nr <- gsub("^.+[[:space:]]", "", head(release)[rel_nr])
  rel_dt <- grep(paste(month.name, collapse = "|"), head(release))
  rel_dt <- gsub("^[[:space:]]+", "", head(release)[rel_dt])
  rel_dt <- mdy(rel_dt)
  cat(silver("Latest release: " %+% magenta$bold(rel_nr) %+% 
               cyan(paste0(" (", rel_dt, ")")) %+% "\n"))
  
  ## Size and name of sequence flat files
  ## ------------------------------------
  file_size <- grep("File Size[[:space:]]+File Name", release)
  id <- grep(paste0("gb", tolower(code), "[[:digit:]]+[.]seq"), release)
  id <- id[id > file_size]
  seq_files <- gsub("^ +", "", release[id])
  seq_files <- do.call(rbind, strsplit(seq_files, "[[:space:]]+"))
  seq_files <- as.data.frame(seq_files, stringsAsFactors = FALSE)
  names(seq_files) <- c("size", "name")
  cat(silver("Number of files: " %+% magenta$bold(nrow(seq_files)) %+% "\n"))
  
  ## check if files have been downloaded and if they are up to date
  ## --------------------------------------------------------------
  path <- file.path(x@params@data.path, "megaptera_data/gb_sequence")
  fn <- file.path(path, paste(seq_files$name, "gz", sep = "."))
  seq_files$local <- sapply(fn, file.exists)
  seq_files$release <- NA
  tmp <- lapply(fn[seq_files$local], gbflatInfo)
  tmp <- data.frame(do.call(rbind, tmp))
  seq_files$release[seq_files$local] <- tmp$release
  cnt1 <- length(which(seq_files$release[seq_files$local] %in% rel_nr))
  cat(silver(" > " %+% magenta$bold(cnt1) %+% " files have a local copy and are up to date\n"))
  cnt2 <- length(which(!seq_files$release[seq_files$local] %in% rel_nr))
  cat(silver(" > " %+% magenta$bold(cnt2) %+% " files have a local copy, but need update\n"))
  cnt3 <- length(which(!seq_files$local))
  cat(silver(" > " %+% magenta$bold(cnt3) %+% " files have no local copy and require download\n"))
  
  ## Download flat files (if any is missing)
  ## ---------------------------------------
  seq_files <- seq_files[!seq_files$local | seq_files$release != rel_nr, ]
  if (!nrow(seq_files)){
    cat(green("All files have a local copy and are up to date\n"))
  } else {
    cat(silver("In total, " %+% magenta$bold(nrow(seq_files)) %+% " files will be downloaded\n"))
    for (i in seq_files$name[1:nrow(seq_files)]){
      download.file(file.path("ftp://ftp.ncbi.nlm.nih.gov/genbank", paste0(i, ".gz")),
                    file.path(path, paste0(i, ".gz")), timeout = 300)
    }
  }
  
  taxa <- dbReadTaxonomy(x)
  taxa <- taxa$taxon[taxa$rank == x@taxon@tip.rank & taxa$status == "scientific name"]

  ## 3. In update modus identify flatfiles that contain the taxon
  ##    that should be updated
  ## -------------------------
  if (!missing(update)){
    taxa <- taxa[grep(update, taxa)]
    cl <- makeCluster(x@params@cpus)
    clusterEvalQ(cl, library(megaptera))
    clusterExport(cl, c("taxa", "fn"), environment())
    id <- parSapply(cl, X = fn, FUN = findInFlatfile, taxa = unlist(taxa)) # 12:13
    stopCluster(cl)
    
    fn <- names(id)[id]
  }
  
  ## 4. Extract sequences from flatfiles and write them to database
  ## --------------------------------------------------------------
  if (length(fn) < x@params@cpus | !x@params@parallel){
    lapply(fn[198:204], gbflat2db, taxa = taxa, megProj = x)
  } else {
    cl <- makeCluster(x@params@cpus)
    clusterEvalQ(cl, library(megaptera))
    clusterExport(cl, c("x", "fn", "taxa"), environment())
    parLapply(cl, X = fn, fun = gbflat2db, taxa = taxa, megProj = x)
    stopCluster(cl)
  }
  
  cat(silver("\nSTEP GenBankFTP finished"))
  td <- Sys.time() - start
  cat(silver(" after " %+% cyan(paste(round(td, 2), attr(td, "units"))) %+% "\n"))
}
