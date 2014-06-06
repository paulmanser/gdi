### function for simulation study of PCA dimension-reduced
### cross covariance test

ccaTest <- function(object, npcs = 3){
  
  if (!is(object, 'GDIset')) stop("object must be a 'GDIset'")
  
  # center and replace NAs with zero for now -------------------------
#   set1.dat <- t(apply(object@set1@dat, 1, na.replace))
#   set2.dat <- t(apply(object@set2@dat, 1, na.replace))
  set1.dat <- object@set1@dat - rowMeans(object@set1@dat)
  set2.dat <- object@set2@dat - rowMeans(object@set2@dat)
  set1.dat[is.na(set1.dat)] <- 0
  set2.dat[is.na(set2.dat)] <- 0
  
  # do PCA grouped by gene -------------------------------------------
  # PM: maybe try to speed this step up
  pc2 <- function(x) prcomp(t(x[, (-1)]))$x[ , 1:min(npcs, nrow(x))]

  # PM: use something else besides dplyr like bigmemory?
  # get PC scores for set 1
  set1.df <- data.frame(gs = object@set1@annot$entrez.id, set1.dat)
  set1.pca <- set1.df %.%
    group_by(gs) %.%
    do2(pc2)
  
  # get PC scores for set 2
  set2.df <- data.frame(gs = object@set2@annot$entrez.id, set2.dat)
  set2.pca <- set2.df %.%
    group_by(gs) %.%
    do2(pc2)
  
  # compute canonical correlation
  cc.results <- mapply(cancor, set1.pca, set2.pca)
  
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


do2 <- function (.data, .f, ...) {
  if (is.null(attr(.data, "indices"))) {
    .data <- dplyr:::grouped_df_impl(.data, attr(.data, "vars"), 
                                     attr(.data, "drop"))
  }
  index <- attr(.data, "indices")
  out <- vector("list", length(index))
  for (i in seq_along(index)) {
    subs <- .data[index[[i]] + 1L, , drop = FALSE]
    out[[i]] <- .f(subs, ...)
  }
  nms <- as.character(attr(.data, "labels")[[1]])
  setNames(out, nms)
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

# replace NAs with mean for now
na.replace <- function(x){
  x.mean <- mean(na.omit(x))
  x[is.na(x)] <- x.mean
  return(x)
}

