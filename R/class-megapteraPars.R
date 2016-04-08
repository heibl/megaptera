## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-02-01)

setClass("megapteraPars", 
         representation = list(
           parallel = "logical",
           cpus = "numeric",
           cluster.type = "character",
           update.seqs = "character", # step B
           retmax = "numeric", # step B
           max.gi.per.spec = "numeric", # step B
           max.bp = "numeric", # step B
           reference.max.dist = "numeric", # step C
           min.seqs.reference = "numeric", # step D
           fract.miss = "numeric", # step G
           filter1 = "numeric",  # step G
           filter2 = "numeric", # step G
           filter3 = "numeric", # step G
           filter4 = "numeric",  # step G
           block.max.dist = "numeric", # step G
           min.n.seq = "numeric", # step G
           gb1 = "numeric", # step F + G
           gb2 = "numeric", # step F + G
           gb3 = "numeric", # step F + G
           gb4 = "numeric", # step F + G
           gb5 = "character") # step F + G
)

"megapteraPars" <- function(...){
  
  params <- list(parallel = FALSE,
                 cpus = 0,
                 cluster.type = "none",
                 update.seqs = "all", 
                 retmax = 500,
                 max.gi.per.spec = 100, 
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
            cat("MEGAPTERA pipeline parameters:",
                "\n             parallel =", object@parallel,   
                "\n                 cpus =", object@cpus,    
                "\n         cluster.type =", object@cluster.type,
                "\n          update.seqs =", object@update.seqs,   
                "\n               retmax =", object@retmax,    
                "\n      max.gi.per.spec =", object@max.gi.per.spec,
                "\n               max.bp =", object@max.bp,             
                "\n   reference.max.dist =", object@reference.max.dist,
                "\n   min.seqs.reference =", object@min.seqs.reference,
                "\n           fract.miss =", object@fract.miss,        
                "\n              filter1 =", object@filter1,             
                "\n              filter2 =", object@filter2,             
                "\n              filter3 =", object@filter3,             
                "\n              filter4 =", object@filter4,            
                "\n       block.max.dist =", object@block.max.dist,     
                "\n            min.n.seq =", object@min.n.seq,           
                "\n                  gb1 =", object@gb1,               
                "\n                  gb2 =", object@gb2,                
                "\n                  gb3 =", object@gb3,               
                "\n                  gb4 =", object@gb4,                 
                "\n                  gb5 =", object@gb5
            )
          }
)