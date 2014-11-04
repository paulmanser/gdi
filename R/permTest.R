

permTest <- function(object, cca.results = NULL, n.perm = 1000, half.life = 300){

  if (!is(object, 'GDIset')) stop("object must be a 'GDIset'")
  
  if (!is.null(cca.results)){
    cat('Performing permutation test on communalities \n')
    incl.ind <- which(unlist(lapply(cca.res$testing.results, function(x) !is.na(x[1]))))
    unique.ids <- names(cca.results$testing.results)[incl.ind]
  }else{
    cat('Performing permutation test on R-squared values \n')
    set1.names <- names(table(getAnnot(object)[[1]]$entrez.id) > 3)
    set2.names <- names(table(getAnnot(object)[[2]]$entrez.id) > 3)
    unique.ids <- intersect(getAnnot(object)[[1]]$entrez.id,
                            getAnnot(object)[[2]]$entrez.id)
  }
    
  out <- foreach(gene=unique.ids, .packages='gdi', .combine='c') %dopar% {
    
    # get locations of sites
    set1.ind <- which(object@set1@annot$entrez.id %in% gene)
    set1.loc <- (end(object@set1@annot[set1.ind]) + start(object@set1@annot[set1.ind]))/2
    names(set1.loc) <- names(object@set1@annot[set1.ind])
    
    set2.ind <- which(object@set2@annot$entrez.id %in% gene)
    set2.loc <- (end(object@set2@annot[set2.ind]) + start(object@set2@annot[set2.ind]))/2
    names(set2.loc) <- names(object@set2@annot[set2.ind])
    
    if (!is.null(cca.results)){
      # save communalities with short variable names
      set1.comm <- cca.results$set1.loadings[[gene]][, 1, drop=FALSE]^2
      set2.comm <- cca.results$set2.loadings[[gene]][, 1, drop=FALSE]^2
      
      # get outer product of communalities
      comm.outer <- set1.comm %*% t(set2.comm)
    
    }else{
      # compute R^2 matrix from actual data
      comm.outer <- cor(t(object@set1@dat[set1.ind, ]),
                        t(object@set2@dat[set2.ind, ]))^2
    }    
    
    # get distance matrix 
    dist.mat <- matrix(set1.loc, nr = length(set1.loc),
                       nc=length(set2.loc))
    
    dist.mat <- abs(sweep(dist.mat, 2, set2.loc, '-'))
    
    # apply exponential decay function to get weights
    lambda <- log(2)/half.life
    weight.mat <- exp(-dist.mat*lambda)
    obs.stat <- sum(comm.outer * weight.mat)
    
    # permute weights and compute stats
    perm.stats <- numeric(n.perm)
    for(jj in 1:n.perm){
      perm.stats[jj] <- sum(comm.outer * weight.mat[sample(nrow(weight.mat)), ][ , sample(ncol(weight.mat))])
    }
    
    perm.pval <- sum(perm.stats > obs.stat)/n.perm
    
    perm.pval
  }
  
  names(out) <- unique.ids
  out  

}

