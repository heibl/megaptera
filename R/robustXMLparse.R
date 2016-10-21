## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-09-15)

robustXMLparse <- function(url, trial = 3, logfile){
  
  ## 1. try to parse XML trial times
  ## -------------------------------
  for (i in 1:trial ){
    err <- try(silent = TRUE,
               xml <- xmlTreeParse(getURL(url), getDTD = FALSE, 
                                   useInternalNodes = TRUE))
    if ( "try-error" %in% class(err) ){
      slog("\n.. trial", i, "FAILED:", err, file = logfile)
    } else break
  } 
  if ( "try-error" %in% class(err) ){
    warning(err)
    return(NULL)
  }
  # saveXML(xml, "robustXMLparse-DEBUG.xml"); system("open -t robustXMLparse-DEBUG.xml")
  
  ## 2. return NULL if 'resources temporarily unavailable'
  ## -----------------------------------------------------
  for ( i in 1:trial ){
    err <- xpathSApply(xml, "/GBSet/Error", xmlValue)
    if ( length(err) > 0 ) {
      slog("\n.. trial", i, "FAILED:", err, file = logfile)
    } else break 
  }
  if ( length(err) > 0 ) {
    warning(err)
    return(NULL)
  }
  
  ## 3. return NULL when results are empty
  ## -------------------------------------
  ERROR <- xpathSApply(xml, fun = xmlValue, path = "//ERROR")
  if ( length(ERROR) > 0 ){
    slog("\n.. FAILED:", ERROR, "..", file = logfile)
    return(NULL)
  } 
  xml
}