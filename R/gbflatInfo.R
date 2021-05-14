#' @export

gbflatInfo <- function(file){
  
  conn <- file(file, open = "r")
  hdr <- readLines(con = conn, n = 10)
  close(conn)
  
  fn <- tolower(gsub("[[:space:]].+$", "", hdr[1]))
  nr <- grep("NCBI-GenBank Flat File Release", hdr)
  nr <- gsub("^.+[[:space:]]", "", head(hdr)[nr])
  dt <- grep(paste(month.name, collapse = "|"), hdr)
  dt <- gsub("^[[:space:]]+", "", head(hdr)[dt])
  dt <- mdy(dt)
  list(file = fn,
       release = nr,
       date = dt)
}

