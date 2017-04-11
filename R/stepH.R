## This code is part of the megaptera package
## © C. Heibl 2016 (last update 2017-02-23)

#' @export
#' @import DBI snowfall

stepH <- function(x, max.mad){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if (!inherits(x, "megapteraProj"))
    stop("'x' is not of class 'megapteraProj'")
  if (x@locus@kind == "undefined") stop("undefined locus not allowed")
  
  ## check if previous step has been run
  ## -----------------------------------
  status <- dbProgress(x)
  if (status$step_g == "pending") {
    stop("the previous step has not been called yet")
  }
  if (status$step_g == "error") {
    stop("the previous step has terminated with an error")
  }
  if (status$step_g == "failure") {
    slog("\nNo data from upstream available - quitting", file = "")
    dbProgress(x, "step_h", "failure")
    return()
  }
  if (status$step_g == "success") {
    dbProgress(x, "step_h", "error")
  }
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tip.rank <- x@taxon@tip.rank
  msa.tab <- paste(tip.rank, gsub("^_", "", gene), sep = "_")
  block.max.dist <- x@params@block.max.dist
  min.n.seq <- x@params@min.n.seq
  max.bp <- x@params@max.bp * 1.5
  if (missing(max.mad)) max.mad <- x@params@max.mad
  
  ## iniate logfile
  ## --------------
  logfile <- paste0("log/", gene, "-stepH.log")
  if (file.exists(logfile)) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste0("\n", Sys.time()),
       "\nSTEP H: detection and separation of unsaturated blocks\n", 
       paste("\n.. locus:", x@locus@sql), file = logfile)
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
  
  ## check if stepC has been run
  ## ---------------------------
  # status <- paste("SELECT DISTINCT status",
  #                 "FROM", acc.tab)
  # status <- dbGetQuery(conn, status)
  # if ("raw" %in% status$status){
  #   dbDisconnect(conn)
  #   stop("stepC has not been run")
  # }
  
  ## check if msa table exists
  ## -------------------------
  # if (!dbExistsTable(conn, msa.tab)){
  #   dbDisconnect(conn)
  #   slog("\nWARNING: table", msa.tab, "does not exist!\n", file = logfile)
  #   td <- Sys.time() - start
  #   slog("\nSTEP H finished after", round(td, 2), attr(td, "units"), "\n",
  #        file = logfile)
  #   return()
  # }
  
  ## check if at least 3 (ingroup) species are available
  ## ---------------------------------------------------
  # n <- paste("SELECT count(spec) FROM", msa.tab)
  # n <- dbGetQuery(conn, n)$count
  # if (n < 3){
  #   dbDisconnect(conn)
  #   slog("\nWARNING: only", n, "species available,", 
  #        "no saturation assessment possible\n", file = logfile)
  #   td <- Sys.time() - start
  #   slog("\nSTEP H finished after", round(td, 2), attr(td, "units"), 
  #        "\n", file = logfile)
  #   dbProgress(x, "step_h", "failure")
  #   return()
  # }
  # if (n < 100){ # 100 is arbitrary
  #   n <- paste("SELECT spec FROM", msa.tab)
  #   n <- dbGetQuery(conn, n)
  #   n <- which(is.ingroup(x, n$spec))
  #   if (length(n) < 3){
  #     dbDisconnect(conn)
  #     slog("\nWARNING:", length(n), "ingroup species available,", 
  #          "no saturation assessment possible\n", file = logfile)
  #     td <- Sys.time() - start
  #     slog("\nSTEP GG finished after", round(td, 2), attr(td, "units"), 
  #          "\n", file = logfile)
  #     dbProgress(x, "step_h", "failure")
  #     return()
  #   }
  # }
  
  ## read alignment
  ## --------------
  slog("\n.. reading alignment with ", file = logfile)
  # if ( x@params@parallel ){
  #   spec <- dbGetQuery(conn, paste("SELECT", tip.rank, "FROM", msa.tab))[, tip.rank]
  #   spec <- paste("^", spec, "$", sep = "")
  #   n <- length(spec)
  #   id <- seq(from = 1, to = n, by = ceiling(n/x@params@cpus))
  #   id <- data.frame(from = id, to = c(id[-1] - 1, n))
  #   spec <- apply(id, 1, function(i, a) paste(a[i[1]:i[2]], collapse = "|"), a = spec)
  #   sfInit(parallel = TRUE, cpus = x@params@cpus, 
  #          type = x@params@cluster.type)
  #   sfLibrary("megaptera", character.only = TRUE)
  #   sfLibrary("seqinr", character.only = TRUE)
  #   megProj <- x
  #   sfExport("spec", "msa.tab", "megProj")
  #   a <- sfLapply(x = megProj, fun = dbPReadDNA, tab.name = msa.tab, 
  #                 regex = TRUE, max.bp = max.bp)
  #   sfStop()
  #   a <- do.call(rbind, a)
  # } else {
  a <- dbReadDNA(x, msa.tab)
  # }
  ## check if stepG has been run properly
  ## should include checking for status = 'raw'
  if (!is.matrix(a)){
    stop("stepG has not been called yet")
  }
  ## how many species are in the table/alignment
  slog(nrow(a), "species", file = logfile)
  
  ## Set status column to 'aligned':
  ## If this is not done and the alignment is below
  ## MAD threshold, blocks from previous runs will
  ## not be eliminated.
  ## -----------------
  SQL <- paste("UPDATE", msa.tab,
               "SET status = 'aligned'",
               "WHERE status ~ 'block'")
  dbSendQuery(conn, SQL)
  
  ## checking saturation of alignment as measures
  ## by the median average difference (MAD) between
  ## corrected and JC69-corrected patristic distances
  ## ------------------------------------------------
  check.mad <- MAD(a)
  slog("\n.. assessing saturation: MAD =", round(check.mad, 5),
       file = logfile)
  if (check.mad <= max.mad){
    SQL <- paste("UPDATE", msa.tab,
                 "SET status = '1 block'",
                 "WHERE status = 'aligned'")
    dbSendQuery(conn, SQL)
    dbDisconnect(conn)
    slog(paste("\n.. MAD below threshold of ", max.mad, ": ",
               "alignment not saturated\n", sep = ""), file = logfile)
    td <- Sys.time() - start
    slog("\nSTEP H finished after", round(td, 2), attr(td, "units"), "\n",
         file = logfile)
    dbProgress(x, "step_h", "success")
    return()
  }
  
  ## prepare guide tree 
  ## ------------------
  slog("\n.. preparing guidetree", file = logfile)
  gt <- comprehensiveGuidetree(x, tip.rank = tip.rank, subset = a)
  gt <- fixNodes(gt)
  a <- a[match(gt$tip.label, rownames(a)), ]
  a <- list(a)
  this.root <- Ntip(gt) + 1
  
  slog("\n.. successively splitting up alignment", file = logfile)
  repeat {
    aa <- splitAlignment(this.root, gt, a[[1]])
    a <- c(a[-1], aa)
    check.mad <- sapply(a, MAD)
    slog("\n   - maximum MAD:", max(check.mad), file = logfile)
    if ( all(check.mad <= max.mad) ) break
    a <- a[order(check.mad, decreasing = TRUE)]
    this.root <- noi(gt, rownames(a[[1]]))
  }
  slog("\n.. alignment was split into", length(a), 
       "unsaturated blocks", file = logfile)
  check.size <- sapply(a, nrow)
  a <- a[order(check.size, decreasing = TRUE)]
  
  ## write alignment blocks to database
  ## ----------------------------------
  tips <- lapply(a, function(z) gsub("_", " ", rownames(z)))
  names(tips) <- paste("block", 1:length(tips))
  tips <- lapply(tips, wrapSQL, term = "taxon", operator = "=", boolean = "OR")
  SQL <- paste("UPDATE", msa.tab,
               "SET", wrapSQL(names(tips), "status", "=", NULL),
               "WHERE", tips)
  lapply(SQL, dbSendQuery, conn = conn)
  
  ## write files
  ## -----------
  #   slog("\n.. write alignment to files ..", file = logfile)
  #   write.phy(a, paste(gene, "masked.phy", sep = "-"))
  #   rownames(a) <- gsub("-", "_", rownames(a))
  #   write.nex(a, paste(gene, "masked.nex", sep = "-"))
  
  ## summary of result
  ## -----------------
  #     slog(paste("\n\n--- final alignment of", gene, "---"),
  #          paste("\nnumber of sequences     :", nrow(a)),
  #          paste("\nnumber of base pairs    :", ncol(a)),
  #          file = logfile)
  #     if ( length(blk) > 1 ) {
  #       slog("\nnumber of blocks        :", length(blk), 
  #            "(", paste(blk, collapse = ", "), ")", 
  #            file = logfile) 
  #     }
  #     slog(paste("\nmean absolute deviation :", this.mad), 
  #          file = logfile)
  #   }
  dbDisconnect(conn)
  slog("\n\nSTEP H finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), "\n", file = logfile)
  dbProgress(x, "step_h", "success")
}
