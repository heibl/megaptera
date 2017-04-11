## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-12-07)

#' @export

taxonomy2html <- function(tax){
  
  if (length(grep("^spec$", names(tax))) == 1){
    ## isolate epitheta
    tax$spec <- gsub("^.*[_| ]", "", tax$spec)
    tip.rank <- "^spec$"
  } else {
    tip.rank <- "^gen$"
  }
  id <- apply(tax[, 1:grep(tip.rank, names(tax))], 2, 
              function(z) length(unique(z)))
  id <- max(which(cumsum(id) == 1:length(id))) -1
  tax <- tax[, -(1:id)]
  
  gen <- which(names(tax) == "gen")
  
  
  cc <- c("<colgroup>",
          paste("<col span=", gen - 1," style='background-color:#ccccff'>", sep = "'"),
          paste("<col span=", 2," style='background-color:#ffff99'>", sep = "'"),
          # paste("<col span=", 2," style='background-color:red'>", sep = "'"),
          "</colgroup>")
  h <- paste("<th>", names(tax), "</th>", sep = "")
  h <- paste(h, collapse = "")
  h <- paste("<tr>", h, "</tr>", sep = "")
  tab <- apply(tax, c(1, 2), function(z) 
    paste("<td>", z, "</td>", sep = ""))
  tab <- apply(tab, 1, paste, collapse = "")
  tab <- paste("<tr>", tab, "</tr>", sep = "")
  tab <- c("<table border = 1>", cc, h, tab, "</table>")
  
  
  z <- c("<!DOCTYPE html>\n\n<html xmlns='http://www.w3.org/1999/xhtml'>", 
         "</head>",
         "<style>", 
         "table, th, td {",
         "border: 1px solid black;",
         " border-collapse: collapse;",
         "}",
         "th, td {",
         "padding: 2px;",
         "}",
         "</style>", "</head>",
         "<body>",
         "<div id='header'>", 
         "<h1 class='title'>Taxonomic classification</h1>",
         tab, "</div>", "</body>", "</html>")
  write(z, file = "report/taxonomy.html")
}