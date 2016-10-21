## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2015-12-18)

dropExtendedIngroup <- function(x, phy, surrogate = FALSE){
  
  extended <- dbReadTaxonomy(x, tag = "extended ingroup", 
                             subset = phy)$spec
  
  if ( surrogate ){
    
    ## delete previous surrogate definitions:
    conn <- dbconnect(x)
    SQL <- paste("UPDATE taxonomy",
                      "SET tag = regexp_replace(tag, '(ingroup [(].+[)])( surrogate:.+$)', '\\1')")
    dbSendQuery(conn, SQL)
    dbDisconnect(conn)
    
    ## case 1: one species of a missing genus
    ingroup <- dbReadTaxonomy(x, tag = "^ingroup")$spec
    present.spec <- intersect(ingroup, phy$tip.label)
    present.gen <- strip.spec(present.spec)
    missing.spec <- setdiff(ingroup, phy$tip.label)
    missing.spec <- data.frame(gen = strip.spec(missing.spec), 
                              spec = missing.spec, 
                              stringsAsFactors = FALSE)
    missing.spec <- split(missing.spec$spec, missing.spec$gen)
    missing.gen <- setdiff(names(missing.spec), present.gen)
    if ( length(missing.gen) > 0 ){
      ## quick and dirty
      id <- grep(paste(missing.gen, collapse = "|"), extended)
      s <- extended[id]
      s <- data.frame(gen = strip.spec(s), spec = s, stringsAsFactors = FALSE)
      s <- split(s$spec, s$gen)
      r <- sapply(missing.spec[names(missing.spec) %in% names(s)], head, n = 1)
      s <- sapply(s, head, n = 1)
      s <- data.frame(spec = sort(r), 
                      surrogate = sort(s), 
                      stringsAsFactors = FALSE)
      apply(s, 1, dbUpdateSurrogate, x = x)
      phy$tip.label[match(s$surrogate, phy$tip.label)] <- s$spec
      extended <- setdiff(extended, s$surrogate)
    }
    
    ## case 2: one species of one congeneric
    present.spec <- intersect(ingroup, phy$tip.label)
    missing.spec <- setdiff(ingroup, phy$tip.label)
    missing.spec <- data.frame(gen = strip.spec(missing.spec), 
                               spec = missing.spec, 
                               stringsAsFactors = FALSE)
    missing.spec <- split(missing.spec$spec, missing.spec$gen)
    one.sister <- sapply(names(missing.spec), 
                         function(p, z) length(grep(p, z)) == 1, 
                         z = present.spec)
    if ( any(one.sister) ){
      ## quick and dirty
      id <- grep(paste(names(one.sister[one.sister]), collapse = "|"), extended)
      s <- extended[id]
      s <- data.frame(gen = strip.spec(s), spec = s, stringsAsFactors = FALSE)
      s <- split(s$spec, s$gen)
      r <- sapply(missing.spec[names(missing.spec) %in% names(s)], head, n = 1)
      s <- sapply(s, head, n = 1)
      s <- data.frame(spec = sort(r), 
                      surrogate = sort(s), 
                      stringsAsFactors = FALSE)
      apply(s, 1, dbUpdateSurrogate, x = x)
      phy$tip.label[match(s$surrogate, phy$tip.label)] <- s$spec
      extended <- setdiff(extended, s$surrogate)
    }
  }
  
  if ( length(extended) == 0 ){
    return(phy)
  }
  phy <- drop.tip(phy, extended)
  phy <- fixNodes(phy)
  phy
}