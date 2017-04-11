## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2017-02-20)

#' @include megapteraPars-class.R
#' @importFrom methods new
#' @export

"megapteraPars" <- function(...){
  
  params <- list(parallel = FALSE,
                 cpus = 0,
                 cluster.type = "none",
                 update.seqs = "all", 
                 retmax = 500,
                 max.gi.per.spec = 10, 
                 max.bp = 5000, 
                 reference.max.dist = 0.25,
                 min.seqs.reference = 10,
                 fract.miss = .25, 
                 filter1 = .5,  # step G
                 filter2 = .25, # step G
                 filter3 = .05, # step G
                 filter4 = .2,  # step G
                 block.max.dist = .5, # step G
                 min.n.seq = 5, # step G
                 max.mad = .01,
                 gb1 = .5, # step F + G
                 gb2 = .5, # step F + G
                 gb3 = 9999, # step F + G
                 gb4 = 2, # step F + G
                 gb5 = "a" # step F + G
  )
  
  args <- list(...)
  notDef <- setdiff(names(args), names(params))
  if ( length(notDef) ) 
    stop ("parameter '", notDef[1], "' is not defined", sep = "")
  
  id <- match(names(args), names(params))
  params[id] <- args
  
  if ( params$parallel & params$cpus == 0 ){
    stop("number of CPUs must be given")
  }
  if ( params$parallel & params$cluster.type == "none"){
    stop("type of cluster must be given: 'SOCK', 'MPI', 'PVM', or 'NWS'")
  }
  
  new("megapteraPars",
      parallel = params$parallel,
      cpus = params$cpus,
      cluster.type = params$cluster.type,
      update.seqs = params$update.seqs, 
      retmax = params$retmax,
      max.gi.per.spec = params$max.gi.per.spec, 
      max.bp = params$max.bp,
      reference.max.dist = params$reference.max.dist,
      min.seqs.reference = params$min.seqs.reference,
      fract.miss = params$fract.miss, 
      filter1 = params$filter1,
      filter2 = params$filter2,
      filter3 = params$filter3,
      filter4 = params$filter4,
      block.max.dist = params$block.max.dist,
      min.n.seq = params$min.n.seq,
      max.mad = params$max.mad,
      gb1 = params$gb1,
      gb2 = params$gb2,
      gb3 = params$gb3,
      gb4 = params$gb4,
      gb5 = params$gb5
  )
}

setMethod("show",
          signature(object = "megapteraPars"),
          function (object) 
          {
            out <- sapply(slotNames(object), slot, object = object)
            names(out) <- format(names(out), justify = "right")
            out <- paste("\n ", names(out), "=", out)
            cat("MEGAPTERA pipeline parameters:", out)
          }
)