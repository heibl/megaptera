## This code is part of the megaptera package
## © C. Heibl 2019 (last update 2019-09-11)

#' @export

megapteraCSS <- function(file){
  
  ## Do not replace an existing CSS file
  ## -----------------------------------
  if (file.exists(file)) return()
  
  ## Define CSS file
  ## ---------------
  css <- c(
    "h1 {",
    "  font-size: 36px;",
    "}",
    "h2 {",
    "  font-size: 24px;",
    "}",
    "table, th, td {",
    "  border: 1px solid black;",
    "  border-collapse: collapse;",
    "}",
    "th, td {",
    "  padding: 8px;",
    "}",
    ".pending {",
    "  width: 12px;", 
    "  height: 12px;",
    "  float: left;",
    "  background: #2a73a4;", 
    "  -moz-border-radius: 6px;", 
    "  -webkit-border-radius: 6px;", 
    "  border-radius: 6px; ",
    "}",
    ".success {", 
    "  width: 12px;", 
    "  height: 12px;", 
    "  float: left;",
    "  background: #007705;", 
    "  -moz-border-radius: 6px;", 
    "  -webkit-border-radius: 6px;", 
    "  border-radius: 6px;",
    "}",
    ".failure {", 
    "  width: 12px;", 
    "  height: 12px;", 
    "  float: left;",
    "  background: #aaaaaa;", 
    "  -moz-border-radius: 6px;", 
    "  -webkit-border-radius: 6px;", 
    "  border-radius: 6px;", 
    "}",
    ".error {",
    "  width: 12px;", 
    "  height: 12px;", 
    "  float: left;",
    "  background: #d30000;", 
    "  -moz-border-radius: 6px;", 
    "  -webkit-border-radius: 6px;", 
    "  border-radius: 6px;", 
    "}"
  )
  write(css, file, sep = "\n")
}






