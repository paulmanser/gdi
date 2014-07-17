
# Define ------------------------------------------------------------------

#' @exportClass GDIset
setClass("GDIset", slots = c(set1 = "GDset", set2 = "GDset"))

# Validate ----------------------------------------------------------------

.validGDIset <- function(object) {
  if (!identical(object@set1@pheno, object@set2@pheno))
    stop("'pheno' must match between GDsets", call. = FALSE)  
  return(TRUE)
}

setValidity("GDIset", .validGDIset)

# Constructors ------------------------------------------------------------

GDIset <- function(x, y)
  new("GDIset", set1 = x, set2 = y)

# Accessors ---------------------------------------------------------------

#' @exportMethod getPheno
setMethod("getPheno", "GDIset", function(object) object@set1@pheno)

#' @exportMethod getDat
setMethod("getDat", "GDIset", function(object){
  out <- list(object@set1@dat, object@set2@dat)
  names(out) <- c(object@set1@platform, object@set2@platform)
  out
})

#' @exportMethod getPlatform
setMethod("getPlatform", "GDIset", function(object){
  list(set1 = object@set1@platform, set2 = object@set2@platform)})

#' @exportMethod getAnnot
setMethod("getAnnot", "GDIset", function(object){
  out <- list(object@set1@annot, object@set2@annot)
  names(out) <- c(object@set1@platform, object@set2@platform)
  out
})

setMethod("[", c("GDIset", "ANY", "ANY"), function(x, i, j, ..., drop = FALSE){
  
  if (!is(i, 'character')){
    stop('Row index must be a vector of entrez ids')
  } else {
    print('Please make sure you are subsetting by entrez id')
  }
  
  GDset1 <- x@set1[which(x@set1@annot$entrez.id %in% i), j]
  GDset2 <- x@set2[which(x@set2@annot$entrez.id %in% i), j]
    
  new("GDIset", set1 = GDset1, set2 = GDset2)
})

setMethod("[", c("GDIset", "missing", "ANY"), function(x, i, j, ..., drop = FALSE){
  
  GDset1 <- x@set1[, j]
  GDset2 <- x@set2[, j]  
  new("GDIset", set1 = GDset1, set2 = GDset2)
})

setMethod("[", c("GDIset", "ANY", "missing"), function(x, i, j, ..., drop = FALSE){
  
  if (!is(i, 'character')){
    stop('Row index must be a vector of entrez ids')
  } else {
    print('Please make sure you are subsetting by entrez id')
  }
  
  GDset1 <- x@set1[which(x@set1@annot$entrez.id %in% i), ]
  GDset2 <- x@set2[which(x@set2@annot$entrez.id %in% i), ]
    
  new("GDIset", set1 = GDset1, set2 = GDset2)
})

getSet <- function(object, whichset = 1){
  if (!is(object, "GDIset"))
    stop("object must be a 'GDIset'")
  
  if (whichset == 1) return(object@set1)
  if (whichset == 2) return(object@set2)
}


# Summaries ---------------------------------------------------------------

setMethod("show", "GDIset", function(object){
  cat("A GDIset object containing \n\n")
  print(object@set1)
  cat("\n\n")
  print(object@set2)
  
})

setMethod("dim", "GDIset", function(x){
  out <- list(c(loci = nrow(x@set1@dat), samples = ncol(x@set1@dat)),
              c(loci = nrow(x@set2@dat), samples = ncol(x@set2@dat)))
  names(out) <- getPlatform(x)
  out
})

# consolidate GDI set -----------------------------------------------------

consolidate <- function(object){
  
  if (!is(object, "GDIset")) stop("object needs to be a 'GDIset'")
  
  genes <- lapply(getAnnot(object), function(x) x$entrez.id)
  has.both <- intersect(genes[[1]], genes[[2]])
  
  set1.include <- which(genes[[1]] %in% has.both)
  set2.include <- which(genes[[2]] %in% has.both)
  
  set1.annot <- object@set1@annot[set1.include]
  set2.annot <- object@set2@annot[set2.include]

  subs1 <- 1:nrow(object@set1@dat) %in% set1.include
  subs2 <- 1:nrow(object@set2@dat) %in% set2.include
  
  set1.dat <- subset(object@set1@dat, subset=subs1)
  set2.dat <- subset(object@set2@dat, subset=subs2)
  
  rownames(set1.dat) <- rownames(object@set1@dat)[set1.include]
  rownames(set2.dat) <- rownames(object@set2@dat)[set2.include]
  
  GDset1 <- GDset(dat = set1.dat,
                  annot = set1.annot,
                  pheno = object@set1@pheno,
                  platform = object@set1@platform)
  
  GDset2 <- GDset(dat = set2.dat,
                  annot = set2.annot,
                  pheno = object@set2@pheno,
                  platform = object@set2@platform)
  
  GDIset(GDset1, GDset2)   
}

