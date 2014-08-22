### function for simulation study of PCA dimension-reduced
### cross covariance test

ccaTest <- function(object, npcs = 5, min.set1=5, min.set2=3){
  
  if (!is(object, 'GDIset')) stop("object must be a 'GDIset'")

  # mean center and combine data sets -------------------------------
  if (is(object@set1@dat, 'ffdf')){
    set1.df <- ffdfrowcenter(X=object@set1@dat)
    set1.df$set <- ff(factor(rep('set1', nrow(set1.df))))  
  } else {
    set1.df <- object@set1@dat - rowMeans(object@set1@dat, na.rm=TRUE)
    set1.df$set <- factor(rep('set1', nrow(set1.df)))
  }
  
  if (is(object@set2@dat, 'ffdf')){
    set2.df <- ffdfrowcenter(X=object@set2@dat)
    set2.df$set <- ff(factor(rep('set2', nrow(set2.df))))
  } else {
    set2.df <- object@set2@dat - rowMeans(object@set2@dat, na.rm=TRUE)
    set2.df$set <- factor(rep('set2', nrow(set2.df)))
  }
  
  if (is(set1.df, 'ffdf') & is(set2.df, 'ffdf')){
    full.set <- ffdfappend(set1.df, set2.df)
  } else {
    full.set <- rbind(as.data.frame(set1.df), as.data.frame(set2.df))
  }

  # do analysis grouped by gene -------------------------------------
  entrez.ids <- c(object@set1@annot$entrez.id, object@set2@annot$entrez.id)
  ids.in.both <- intersect(object@set1@annot$entrez.id, object@set2@annot$entrez.id)
  unique.ids <- unique(ids.in.both)
  ind <- 1
  out.return <- list()
  
#   out <- foreach(gene = unique.ids, .packages='gdi')  %do% {
  for(gene in unique.ids){
    
    cat(gene, ' ', ind)
    
    dat <- full.set[which(entrez.ids == gene), ]
    n.sites <- table(dat$set)
    
    if (n.sites[1] < min.set1 | n.sites[2] < min.set2){
      out.return[[gene]] <- NA
      cat(' omitted')
    } else {
      #mean center and replace NA's with zero (centered mean) for now
      dat[is.na(dat)] <- 0      
        
      # perform PCA for each set
      pc.out <- (dat %.%
                 group_by(set) %.%
                 do(pcs = pca(.)))$pcs
      
      pcs.1 <- t(pc.out[[1]]$x[, 1:min(npcs, ncol(pc.out[[1]]$x)), drop=FALSE])
      pcs.2 <- t(pc.out[[2]]$x[, 1:min(npcs, ncol(pc.out[[2]]$x)), drop=FALSE])
      
      # do CCA on PC scores
      cc.res <- cancor(t(pcs.1), t(pcs.2))
      
      # do bartlett's significance test for CCA
      n <- ncol(dat) - 1
      npcs1 <- nrow(pcs.1)
      npcs2 <- nrow(pcs.2)
      test.stat <- -(n - 1 - .5*(npcs1+npcs2+1)) * sum(log(1-cc.res$cor^2))
      df <- npcs1 * npcs2
      p.value <- 1 - pchisq(test.stat, df)
      
      # compute redundancy for first CC
      set1.scores <- t(pcs.1) %*% cc.res$xcoef[, 1, drop=FALSE]
      set2.scores <- t(pcs.2) %*% cc.res$ycoef[, 1, drop=FALSE]
      
      set1.load <- cor(set1.scores, t(subset(dat, subset=set=='set1')[, -ncol(dat)]))
      set2.load <- cor(set2.scores, t(subset(dat, subset=set=='set2')[, -ncol(dat)]))
      
      set1.redundancy <- mean(cor(set2.scores, t(subset(dat, subset=set=='set1')[, -ncol(dat)]))^2)
      set2.redundancy <- mean(cor(set1.scores, t(subset(dat, subset=set=='set2')[, -ncol(dat)]))^2)
      
      output <- list()
      output$test.results <- c(chisq_stat=test.stat, df=df, 
                               p_value=p.value, 
                               set1_r2=set1.redundancy, 
                               set2_r2=set2.redundancy)
      
      output$loadings <- list()
      output$loadings$set1 <- set1.load
      output$loadings$set2 <- set2.load
      
      output$scores <- list()
      output$scores$set1 <- set1.scores
      output$scores$set2 <- set2.scores
      
      out.return[[gene]] <- output
    }
    ind <- ind + 1
    cat('\n')
  }

  # consolidate results  ---------------------------------------
  out.final <- list()

  # testing results
  out.final$testing.results <- t(sapply(out.return, function(x){
    if(!is.na(x)){
      x$test.results
    } else {
      rep(NA, 5)
    }
  }))
  
  # set 1 cc scores
  out.final$set1.scores <- t(sapply(out.return, function(x, n.c){
    if(!is.na(x)){
      x$scores$set1
    } else {
      rep(NA, n.c)
    }
  }, n.c=dim(full.set)[2]))

  # set 2 cc scores
  out.final$set2.scores <- t(sapply(out.return, function(x, n.c){
    if(!is.na(x)){
      x$scores$set2
    } else {
      rep(NA, n.c)
    }
  }, n.c=dim(full.set)[2]))

  # set1 loadings
  out.final$set1.loadings <- lapply(out.return, function(x){
    if(!is.na(x)){
      x$loadings$set1
    } else {
      NA
    }
  })

  # set 2 loadings
  out.final$set2.loadings <- lapply(out.return, function(x){
    if(!is.na(x)){
      x$loadings$set2
    } else {
      NA
    }
  })
  out.final

}

ffdfrowcenter <- function(X){    

  stopifnot(is.ffdf(X))    
  xchunks <- chunk(X)    
  result <- NULL    
       
  for (i in xchunks){                     
    res.chunk <- t(as.matrix(X[i, ]) - rowMeans(as.matrix(X[i, ])), na.rm=TRUE)
    res.chunk[is.na(res.chunk)] <- 0
    
    if (!is.null(dim(res.chunk))){
      res.chunk <- t(res.chunk)
      colnames(res.chunk) <- colnames(X)
    }
    res.chunk <- as.ffdf(data.frame(res.chunk))
    result <- ffdfappend(result, res.chunk)
  }

  result
}

pca <- function(x) prcomp(t(x[, -ncol(x)]))
