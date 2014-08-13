### function for simulation study of PCA dimension-reduced
### cross covariance test

ccaTest <- function(object, npcs = 3){
  
  if (!is(object, 'GDIset')) stop("object must be a 'GDIset'")
  
  # center and replace NAs with zero for now -------------------------
<<<<<<< HEAD
  set1.df <- object@set1@dat - rowMeans(object@set1@dat)
  set2.df <- object@set2@dat - rowMeans(object@set2@dat)
  set1.df[is.na(set1.df)] <- 0
  set2.df[is.na(set2.df)] <- 0
  
  # do PCA grouped by gene -------------------------------------------
  # get PC scores and redundancies
  set1.pca <- mclapply(split(set1.df, object@set1@annot$entrez.id), pc2, npcs=npcs)
  set2.pca <- mclapply(split(set2.df, object@set2@annot$entrez.id), pc2, npcs=npcs)
    
=======
  mean.center <- function(x){
    z <- x - mean(x)
    z[is.na(z)] <- 0
    z
  }
  
  set1.df <- ffdfrowapply(X=object@set1@dat, FUN=mean.center)  
  set2.df <- ffdfrowapply(X=object@set2@dat, FUN=mean.center)
  
  set1.df$entrez.id <- ff(factor(object@set1@annot$entrez.id))
  set2.df$entrez.id <- ff(factor(object@set2@annot$entrez.id))
  
  # do PCA grouped by gene -------------------------------------------
  pc2 <- function(x){
    z <- prcomp(t(x[, (-ncol(x))]))$x[ , 1:min(npcs, nrow(x))]
    as.data.frame(z)
  }
  
  pc3 <- function(x) as.data.frame(lapplyBy(~entrez.id, data=x, pc2))
  
  # get PC scores for set 1
  set1.pca <- ffdfdply(x=set1.df, split=set1.df$entrez.id,
                       BATCHBYTES=5e5, FUN=pc3, trace=FALSE)
  
  # get PC scores for set 2
  set2.pca <- ffdfdply(x=set2.df, split=set2.df$entrez.id,
                       BATCHBYTES=5e5, FUN=pc3, trace=FALSE)

  # split PC scores by entrez id
  get.entrez <- function(x) strsplit(x, '\\.')[[1]][1]
  
  set1.entrez <- sapply(colnames(set1.pca), get.entrez)
  set1.pca.list <- lapply(as.list(unique(set1.entrez)), 
                          FUN = function(x, grp, data) data[,grp %in% x],
                          grp = set1.entrez, data = as.data.frame(set1.pca))
  names(set1.pca.list) <- unique(set1.entrez)
  
  set2.entrez <- sapply(colnames(set2.pca), get.entrez)
  set2.pca.list <- lapply(as.list(unique(set2.entrez)),
                          FUN = function(x, grp, data) data[, grp %in% x],
                          grp = set2.entrez, data = as.data.frame(set2.pca))
  names(set2.pca.list) <- unique(set2.entrez)
  
>>>>>>> use_ff
  # compute canonical correlation
  cc.results <- mcmapply(cancor, set1.pca.list, set2.pca.list)
  
  # compute significance test
  test.results <- apply(cc.results, 2, cc.sig.test, n = nrow(getPheno(object)))
<<<<<<< HEAD
  
  # compute redundancy coefs
  redundancy.results <- mapply(cc.redundancy, 
                               cc.res = as.list(as.data.frame(cc.results)),
                               pca.res1 = set1.pca,
                               pca.res2 = set2.pca, 
                               set1.df = split(set1.df, object@set1@annot$entrez.id),
                               set2.df = split(set2.df, object@set2@annot$entrez.id))  
  
  out <- mapply(cbind, test.results, redundancy.results)
  out <- apply(out, 2, as.data.frame)  
  out
}

pc2 <- function(z, npcs) prcomp(t(z))$x[ , 1:min(npcs, nrow(z))]                            
                            
=======

  # compute redundancy coefs

  set1.v <- ffdfrowapply(subset(set1.df, select=1:(ncol(set1.df)-1)), var)
  set1.v <- as.data.frame(set1.v)
  set2.v <- ffdfrowapply(subset(set2.df, select=1:(ncol(set2.df)-1)), var)
  set2.v <- as.data.frame(set2.v)
                                      
  redundancy <- foreach(e.id=colnames(cc.results)) %do% {
    
    set1.scores <- as.matrix(set1.pca.list[[e.id]]) %*% cc.results[, e.id]$xcoef
    set2.scores <- as.matrix(set2.pca.list[[e.id]]) %*% cc.results[, e.id]$ycoef
    
    set1.subs <- which(object@set1@annot$entrez.id %in% e.id)
    set2.subs <- which(object@set2@annot$entrez.id %in% e.id)
    
    set1.dat <- set1.df[set1.subs, 1:(ncol(set1.df)-1)]
    set2.dat <- set2.df[set2.subs, 1:(ncol(set2.df)-1)]
    
    set1.comm  <- cor(t(set1.dat), set1.scores)^2
    set2.comm  <- cor(t(set2.dat), set2.scores)^2
    n.ccs <- length(cc.results[, e.id]$cor)
    
    set1.r2 <- colSums((set1.v[set1.subs] * set1.comm)/sum(set1.v[set1.subs]))
    set2.r2 <- colSums((set2.v[set2.subs] * set2.comm)/sum(set2.v[set2.subs]))
    set1.redundancy <- set1.r2[1:n.ccs]*(cc.results[, e.id]$cor^2)
    set2.redundancy <- set2.r2[1:n.ccs]*(cc.results[, e.id]$cor^2)
    out <- cbind(set1.redundancy, set2.redundancy)
    out    
  }

  out <- mapply(cbind, test.results, redundancy)
  out <- apply(out, 2, as.data.frame)  
  out
}


>>>>>>> use_ff
# lrt for each of CC pairs and choose how many CCs to keep
cc.sig.test <- function(object, n){
  npcs.set1 <- nrow(object[2]$xcoef)
  npcs.set2 <- nrow(object[3]$ycoef)
  cc.rho2 <- rev(object[1]$cor^2)
  chisq.stat <- rev((-1) * (n - 1 - .5 * (npcs.set1 + npcs.set2 + 1)) * log(cumprod( 1 - cc.rho2 ) ))
  df <- (npcs.set1 - 1:length(cc.rho2) + 1 ) * (npcs.set2 - 1:length(cc.rho2) + 1 )
  p.values <- (1 - pchisq(chisq.stat, df))
  results.df <- data.frame(chisq.stat, df, p.values)  
}



# compute redundancy coefs
cc.redundancy <- function(cc.res, pca.res1, pca.res2, set1.df, set2.df){
  # add redundancy coefficient  
  set1.scores <- pca.res1 %*% cc.res$xcoef
  set2.scores <- pca.res2 %*% cc.res$ycoef
  set1.comm  <- cor(t(set1.df), set1.scores)^2
  set2.comm  <- cor(t(set2.df), set2.scores)^2
  
  # redundancy index  
  set1.v <- colVars(t(set1.df))
  set2.v <- colVars(t(set2.df))
  n.ccs <- length(cc.res$cor)
  set1.r2 <- colSums((set1.v * set1.comm)/sum(set1.v))
  set2.r2 <- colSums((set2.v * set2.comm)/sum(set2.v))
  set1.redundancy <- set1.r2[1:n.ccs]*(cc.res$cor^2)
  set2.redundancy <- set2.r2[1:n.ccs]*(cc.res$cor^2)
  out <- cbind(set1.redundancy, set2.redundancy)
  return(out)
}

<<<<<<< HEAD
=======
ffdfrowapply <- function(X, FUN){    

  stopifnot(is.ffdf(X))    
  xchunks <- chunk(X)    
  result <- NULL    
       
  for (i in xchunks){                     
    res.chunk <- apply(as.matrix(X[i, ]), 1, FUN)
    if (!is.null(dim(res.chunk))){
      res.chunk <- t(res.chunk)
      colnames(res.chunk) <- colnames(X)
    }
    res.chunk <- as.ffdf(data.frame(res.chunk))
    result <- ffdfappend(result, res.chunk)
  }

  result
}

  
  
  
  
  
>>>>>>> use_ff
