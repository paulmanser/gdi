### function for simulation study of PCA dimension-reduced
### cross covariance test

ccaTest <- function(object, npcs = 3){
  
  if (!is(object, 'GDIset')) stop("object must be a 'GDIset'")
  
  # center and replace NAs with zero for now -------------------------
  set1.df <- object@set1@dat - rowMeans(object@set1@dat)
  set2.df <- object@set2@dat - rowMeans(object@set2@dat)
  set1.df[is.na(set1.df)] <- 0
  set2.df[is.na(set2.df)] <- 0
  
  # do PCA grouped by gene -------------------------------------------
  # get PC scores and redundancies
  set1.pca <- mclapply(split(set1.df, object@set1@annot$entrez.id), pc2, npcs=npcs)
  set2.pca <- mclapply(split(set2.df, object@set2@annot$entrez.id), pc2, npcs=npcs)
    
  # compute canonical correlation
  cc.results <- mapply(cancor, set1.pca, set2.pca)
  
  # compute significance test
  test.results <- apply(cc.results, 2, cc.sig.test, n = nrow(getPheno(object)))
  
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

