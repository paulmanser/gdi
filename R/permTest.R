

permTest <- function(object, min.set1 = 5, min.set2 = 3, n.perm = 1000, lambda = 300){

  if (!is(object, 'GDIset')) stop("object must be a 'GDIset'")
  
  # replace NAs with mean for now
  set1.df <- object@set1@dat - rowMeans(object@set1@dat)
  set2.df <- object@set2@dat - rowMeans(object@set2@dat)
  set1.df[is.na(set1.df)] <- 0
  set2.df[is.na(set2.df)] <- 0
  
  # create data.frames for split/apply
  set1.df <- data.frame(start = start(object@set1@annot),
                        end = end(object@set1@annot),
                        entrez.id = factor(object@set1@annot$entrez.id), 
                        type = factor(rep('set1', nrow(set1.df))),
                        set1.df)
  
  set2.df <- data.frame(start = start(object@set2@annot),
                        end = end(object@set2@annot),
                        entrez.id = factor(object@set2@annot$entrez.id), 
                        type = factor(rep('set2', nrow(set2.df))),
                        set2.df)
  
  int.df <- rbind(set1.df, set2.df) 
  
  # get R^2 matrices 
  int.cov2 <- mclapply(split(int.df, int.df$entrez.id), get.int.cov2)
  
  # get distance matrix using annot
  int.dist <- mclapply(split(int.df, int.df$entrez.id), get.int.dist, lambda=lambda)
  
  # permute!
  perm.list <- mapply(int.cov2, int.dist, FUN=list, SIMPLIFY=FALSE)
  mclapply(perm.list, permute.weights, n.perm = n.perm)  

}


get.int.cov2 <- function(x){
  
  cov(t(as.matrix(x[x$type == 'set1', -(1:4)])),
      t(as.matrix(x[x$type == 'set2', -(1:4)])))^2
  
}

get.int.dist <- function(x, lambda){
  
  # get distance from middle of exon/region for now
  location <- (x$start + x$end)/2
  
  # get distances
  dist.mat <- matrix(location[x$type == 'set1'], 
                     nr = sum(x$type == 'set1'),
                     nc = sum(x$type == 'set2'))
  
  dist.mat <- abs(sweep(dist.mat, 2, location[x$type == 'set2'], '-'))
  
  # apply weight function
  exp(-dist.mat/lambda)
  
}

permute.weights <- function(r2.dist, n.perm, min.set1 = 5, min.set2 = 3){
  
  if (nrow(r2.dist[[1]]) < min.set1  | ncol(r2.dist[[1]]) < min.set2)
    return(NA)
    
  r2 <- r2.dist[[1]]
  weights <- r2.dist[[2]]
  
  perm.stats <- numeric(n.perm)
  orig.stat <- sum(r2 * weights)
  
  # permute rows and columns
  for(jj in 1:n.perm){
    perm.stats[jj] <- sum(r2 * weights[sample(nrow(weights)), ][ , sample(ncol(weights))])
  }
  
  perm.pval <- sum(perm.stats > orig.stat)/n.perm
}



