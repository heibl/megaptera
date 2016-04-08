## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2016-02-15)

htmlTable <- function(tab) {
  
  ## table header
  h <- paste("<th>", colnames(tab), "</th>", sep = "")
  h <- paste(h, collapse = "")
  h <- paste("<tr>", h, "</tr>", sep = "")
  
  ## table body
  tab <- apply(tab, c(1, 2), function(z) 
    paste("<td>", z, "</td>", sep = ""))
  tab <- apply(tab, 1, paste, collapse = "")
  tab <- paste("<tr>", tab, "</tr>", sep = "")
  tab <- c("<table border = 1>", h, tab, "</table>")
  tab
}