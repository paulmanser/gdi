### function for simulation study of PCA dimension-reduced
### cross covariance test

ccaTest <- function(object, npcs = 3){
  
  if (!is(object, 'GDIset')) stop("object must be a 'GDIset'")
  
  # center and replace NAs with zero for now -------------------------
  mean.center <- function(x){
    z <- scale(x, scale=FALSE)
    z[is.na(z)] <- 0
    z
  }

  # need to make ffappply return ff object
  set1.dat <- t(ffapply(X=object@set1@dat, MARGIN=1, AFUN="mean.center", 
                        CFUN='ccbind', RETURN=TRUE))
  entrez.id1 <- as.character(object@set1@annot$entrez.id)
  set1.dat <- as.ffdf(data.frame(set1.dat, entrez.id=entrez.id1))
  
  set2.dat <- t(ffapply(X=object@set2@dat, MARGIN=1, AFUN="mean.center", 
                      RETURN=TRUE, CFUN='crbind'))
  entrez.id2 <- as.character(object@set2@annot$entrez.id)
  set2.dat <- as.ffdf(data.frame(set2.dat, entrez.id=entrez.id2))
  
  # do PCA grouped by gene -------------------------------------------
  # PM: maybe try to speed this step up
  pc2 <- function(x){
    z <- prcomp(t(x[, (-ncol(x))]))$x[ , 1:min(npcs, nrow(x))]
    as.data.frame(z)
  }
  
  pc3 <- function(x) as.data.frame(lapplyBy(~entrez.id, data=x, pc2))
  
  # get PC scores for set 1
  set1.pca <- ffdfdply(x=set1.dat, 
                       split=set1.dat$entrez.id,
                       BATCHBYTES=5e5,
                       FUN=pc3, trace=FALSE)
  
  # get PC scores for set 2
  set2.pca <- ffdfdply(x=set2.dat,
                       split=set2.dat$entrez.id,
                       BATCHBYTES=5e5,
                       FUN=pc3, trace=FALSE)

  get.entrez <- function(x) strsplit(x, '\\.')[[1]][1]
  
  set1.pca.list <- as.data.frame(t(as.data.frame(set1.pca)))
  set1.pca.list <- split(set1.pca.list, 
                         sapply(rownames(set1.pca.list), get.entrez))
    
  set2.pca.list <- as.data.frame(t(as.data.frame(set1.pca)))
  set2.pca.list <- split(set2.pca.list,
                         sapply(rownames(set2.pca.list), get.entrez))  
    
  # compute canonical correlation
  # use cancor and something else here
  cc.results <- mapply(cancor, set1.pca.list, set2.pca.list)
  
  # compute significance test
  test.results <- apply(cc.results, 2, cc.sig.test, n = nrow(getPheno(object)))
  
  # compute redundancy coefs
  redundancy.results <- mapply(cc.redundancy, 
                               cc.res = as.list(as.data.frame(cc.results)),
                               pca.res1 = set1.pca,
                               pca.res2 = set2.pca, 
                               set1.dat = split(set1.dat, object@set1@annot$entrez.id),
                               set2.dat = split(set2.dat, object@set2@annot$entrez.id))  
  
  out <- mapply(cbind, test.results, redundancy.results)
  out <- apply(out, 2, as.data.frame)  
  out
}

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
cc.redundancy <- function(cc.res, pca.res1, pca.res2, set1.dat, set2.dat){
  # add redundancy coefficient  
  set1.scores <- pca.res1 %*% cc.res$xcoef
  set2.scores <- pca.res2 %*% cc.res$ycoef
  set1.comm  <- cor(t(set1.dat), set1.scores)^2
  set2.comm  <- cor(t(set2.dat), set2.scores)^2
  
  # redundancy index  
  set1.v <- colVars(t(set1.dat))
  set2.v <- colVars(t(set2.dat))
  n.ccs <- length(cc.res$cor)
  set1.r2 <- colSums((set1.v * set1.comm)/sum(set1.v))
  set2.r2 <- colSums((set2.v * set2.comm)/sum(set2.v))
  set1.redundancy <- set1.r2[1:n.ccs]*(cc.res$cor^2)
  set2.redundancy <- set2.r2[1:n.ccs]*(cc.res$cor^2)
  out <- cbind(set1.redundancy, set2.redundancy)
  return(out)
}

