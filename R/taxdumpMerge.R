


taxdumpMerge <- function(a, b, linconf = "b"){
  
  
  
  for (i in 3:nrow(b)){
    
    save(a, b, i, file = "user_data/DEV_taxdumpMerge2.rda")
    
    bb <- b[i, ]
    aa <- a[a$taxon == bb$taxon & a$rank == bb$rank, ]
    
    cat(silver(bb$taxon %+% " ... "))
    ## Is taxon name of same rank present?
    if (!nrow(aa)){
      
      ## Find anchor point
      ## -----------------
      bbLIN <- taxdumpLineage(b, bb$taxon)
      aaLINbb <- intersect(bbLIN$taxon, a$taxon)
      aaLINbb <- aaLINbb[aaLINbb != "root"]
      
      anchor <- aaLINbb[which.min(match(aaLINbb, bbLIN$taxon))]
      a <- taxdumpAddNode(a, rank = bb$rank, taxon = bb$taxon, parent = anchor)
      
      cat(red("inserted at node '" %+% bold(anchor) %+% "'\n"))
      next
    }
    
    ## Do the have the same status?
    ## ----------------------------
    if (aa$status != bb$status){
      aaa <- a[a$id == aa$id, ]
      bbb <- b[b$id == bb$id, ]
      aaa <- rbind(aaa, a[a$taxon %in% bbb$taxon, ])
      bbb <- rbind(bbb, b[b$taxon %in% aaa$taxon, ])
      aaa <- unique(rbind(aaa, a[a$id %in% aaa$id, ]))
      bbb <- unique(rbind(bbb, b[b$id %in% bbb$id, ]))
      if (length(setdiff(aaa$taxon, bbb$taxon)) | length(setdiff(bbb$taxon, aaa$taxon))){
        stop("implement me!")
      }
      names(aaa) <- paste(names(aaa), "a", sep = "_")
      names(bbb) <- paste(names(bbb), "b", sep = "_")
      if (length(setdiff(aaa$taxon_a, bbb$taxon_b)))
      ab <- cbind(aaa, bbb[match(aaa$taxon, bbb$taxon), ])
      ## Resolve in favor of 'b'
      if (ab$status_a[ab$taxon_a == bb$taxon] == "synonym" & 
          ab$status_b[ab$taxon_a == bb$taxon] == "scientific name"){
        
        # Set status
        a$status[a$taxon == bb$taxon] <- "scientific name"
        # Care for synonyms
        syns <- ab$taxon_a[ab$id_a == ab$id_a[ab$taxon_a == bb$taxon]]
        syns <- setdiff(syns, bb$taxon)
        if (length(syns) > 1) stop("implement me! [code3]")
        if (length(syns)){
          if (ab$status_b[ab$taxon_b == syns] == "scientific name"){
            ## Synonym gets own node and status "scientific name"
            if (nrow(a[a$parent_id == a$id[a$taxon == syns], ])){
              stop("implement me [split node has children]")
            }
            a$id[a$taxon == syns] <- max(c(a$id, a$parent_id)) + 1
            a$status[a$taxon == syns] <- "scientific name"
          }
        }
        
        
      } else {
        stop("implement me!")
      }
      
    }
    
    ## Are lineages congruent?
    ## They need to have at least one taxon in common and each lineage must not
    ## contain elements that are not present in the other lineage, but in the
    ## other parent-child taxonomy table
    ## ---------------------------------
    aaLIN <- taxdumpLineage(a, bb$taxon)
    aaLIN <- aaLIN[!aaLIN$taxon %in% c("root", bb$taxon), ]
    bbLIN <- taxdumpLineage(b, bb$taxon)
    bbLIN <- bbLIN[!bbLIN$taxon %in% c("root", bb$taxon), ]
    
    aaLINbb <- intersect(aaLIN$taxon, bbLIN$taxon)
    if (!length(aaLINbb)){
      stop("Lineages do not have a common taxon")
    }
    aaPRIVLIN <- setdiff(aaLIN$taxon, bbLIN$taxon)
    if (length(aaPRIVLIN)){
      test <- intersect(aaPRIVLIN, b$taxon)
      if (length(test)){
        stop("Lineages are not congruent: ",  test)
      }
    }
    bbPRIVLIN <- setdiff(bbLIN$taxon, aaLIN$taxon)
    if (length(bbPRIVLIN)){
      test <- intersect(bbPRIVLIN, a$taxon)
      if (length(test)){
        stop("Lineages are not congruent: ", test)
      }
    }
    cat(green("OK\n"))
    
    
    
    
    
    
    
  }
}
