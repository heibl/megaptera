parseTaxonomy <- function(){
  
  NCBI_names <- read.table("names.dmp", sep = "|", 
                      # nrows = 9773,
                      quote = "",
                      strip.white = TRUE, comment.char = "")
  NCBI_names <- NCBI_names[, 1:4]
  names(NCBI_names) <- c("tax_id", "name_txt", "name_unique", "name_class")
  
  conn <- dbConnect(PostgreSQL(), dbname = "taxonomy",
                    user = "postgres", password = "oxalis")
  dbRemoveTable(conn, "ncbi_names")
  SQL <- "CREATE TABLE ncbi_names (
  tax_id integer NOT NULL,
  name_txt text NOT NULL,
  name_unique text NOT NULL,
  name_class text,
  CONSTRAINT ncbi_names_pk PRIMARY KEY (tax_id, name_txt, name_unique)
  )
  WITH (OIDS=FALSE);
  ALTER TABLE ncbi_names OWNER TO postgres;"
  dbSendQuery(conn, SQL)
  dbWriteTable(conn, "ncbi_names", NCBI_names, 
               row.names = FALSE, append = TRUE)
  
  dbDisconnect(conn)
  
}