### function for simulation study of PCA dimension-reduced
### cross covariance test

ccaTest <- function(object, npcs = 5, min.set1=5, min.set2=3, cc.pvalue.threshold=.1){
  
  if (!is(object, 'GDIset')) stop("object must be a 'GDIset'")

  # combine data sets -------------------------------
  set1.df <- object@set1@dat
  if (is(object@set1@dat, 'ffdf')){
    set1.df$set <- ff(factor(rep('set1', nrow(set1.df))))  
  } else {
    set1.df$set <- factor(rep('set1', nrow(set1.df)))
  }
  
  set2.df <- object@set2@dat
  if (is(object@set2@dat, 'ffdf')){
    set2.df$set <- ff(factor(rep('set2', nrow(set2.df))))
  } else {
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
    set1 <- as.matrix(dat[dat$set == 'set1', -ncol(dat)])
    set2 <- as.matrix(dat[dat$set == 'set2', -ncol(dat)])
    
    set1 <- apply(set1, 1, na2mean)
    set2 <- apply(set2, 1, na2mean)
    
    if (n.sites[1] < min.set1 | n.sites[2] < min.set2){
      out.return[[gene]] <- NA
      cat(' omitted')
    } else {
        
      # perform PCA for each set
      pca.set1 <- prcomp(set1)
      pca.set2 <- prcomp(set2)
      
      # do CCA on PC scores
      cc.res <- cancor(pca.set1$x[, 1:npcs], pca.set2$x[, 1:npcs])
      
      # do LRT w/ bartlett correction for CCA
      n <- nrow(set1)
      cc.rho2 <- rev(cc.res$cor^2)
      test.stat <- (-1)*(n - 1 - .5 * (npcs + npcs + 1)) * log(cumprod(1 - cc.rho2))
      df <- (npcs - length(cc.rho2):1 + 1 ) * (npcs - length(cc.rho2):1 + 1 )
      p.value <- (1 - pchisq(test.stat, df))
      
      test.stat <- rev(test.stat)
      df <- rev(df)
      p.value <- rev(p.value)      
      
      n.ccs <- npcs
      
      # compute canonical covariate scores and loadings --------------------------------------
      
      set1.scores <- pca.set1$x[, 1:npcs, drop=FALSE] %*% cc.res$xcoef[, 1:n.ccs, drop=FALSE]
      set2.scores <- pca.set2$x[, 1:npcs, drop=FALSE] %*% cc.res$ycoef[, 1:n.ccs, drop=FALSE]
      set1.loads <- cor(set1, set1.scores)
      set2.loads <- cor(set2, set2.scores)
      
      # redundancy index  
      dat2cc.set1 <- colSums((colVars(set1) * set1.loads^2)/sum(colVars(set1)))
      dat2cc.set2 <- colSums((colVars(set2) * set2.loads^2)/sum(colVars(set2)))
      set1.redundancy <- dat2cc.set1*(cc.res$cor^2)[1:n.ccs]
      set2.redundancy <- dat2cc.set2*(cc.res$cor^2)[1:n.ccs]
      
      # consolidate results into a list --------------------------------------      
      output <- list()
      output$test.results <- cbind(chisq_stat=test.stat, 
                                   df=df, 
                                   p_value=p.value, 
                                   set1_r2=set1.redundancy, 
                                   set2_r2=set2.redundancy)
      
      output$loadings <- list()
      output$loadings$set1 <- set1.loads
      output$loadings$set2 <- set2.loads
      
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
  out.final$testing.results <- lapply(out.return, function(x){
    if(!is.na(x)){
      x$test.results
    } else {
      rep(NA, 5)
    }
  })
  
  # set 1 cc scores
  out.final$set1.scores <- lapply(out.return, function(x, n.c){
    if(!is.na(x)){
      x$scores$set1
    } else {
      rep(NA, n.c)
    }
  }, n.c=dim(full.set)[2])

  # set 2 cc scores
  out.final$set2.scores <- lapply(out.return, function(x, n.c){
    if(!is.na(x)){
      x$scores$set2
    } else {
      rep(NA, n.c)
    }
  }, n.c=dim(full.set)[2])

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

na2mean <- function(x){
  x[is.na(x)] <- mean(na.omit(x))
  return(x)
}
pca <- function(x) prcomp(t(x[, -ncol(x)]))

