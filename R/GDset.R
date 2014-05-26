
# Define ------------------------------------------------------------------

#' @exportClass GDset 
setClass("GDset",
         slots = c(dat = 'data.frame',
                   annot = "GRanges",
                   pheno = 'data.frame',
                   platform = "character"
                   ))

# Validate ----------------------------------------------------------------

.validGDset <- function(object) {
  
  # Required annotation columns
  annot.req <- c(entrez.id = "character")
  
  # Check for required annotation column(s)
  md <- mcols(object@annot)
  if (!all(names(annot.req) %in% names(md))) {
    stop("GDset slot 'annot' must contain all of the following columns:\n", 
         paste(names(columns.req), collapse = "\n"), call. = FALSE)
  }
  
  # Check that annotation matches data
  if (!identical(rownames(object@dat), names(object@annot))){
    stop("Names of 'annot' must match rownames of 'experimentData'", call. = FALSE)
  }
  
  # Check that meta data matches data
  if (!identical(colnames(object@dat), rownames(object@pheno))){
    stop("colnames of 'dat' must match rownames of 'pheno'")
  }
    
  return(TRUE)
}

setValidity("GDset", .validGDset)

# Constructors ------------------------------------------------------------

GDset <- function(dat, annot, pheno, platform){
  new("GDset", 
      dat = dat,
      pheno = pheno,
      platform = platform,
      annot = annot)  
}

# Accessors ---------------------------------------------------------------

#' @export getPlatform
setMethod("getPlatform", "GDset", function(object) object@platform)

#' @export getAnnot
setMethod("getAnnot", "GDset", function(object) object@annot)

#' @export getPheno
setMethod("getPheno", "GDset", function(object) object@pheno)

#' @export getDat
setMethod("getDat", "GDset", function(object) object@dat)

setMethod("[", c("GDset", "ANY", "ANY"),
          function(x, i, j, ..., drop = FALSE){
            new("GDset", annot = x@annot[i], 
                dat = x@dat[i, j, drop=FALSE],
                pheno = x@pheno[j, , drop=FALSE], 
                platform = x@platform)  
          })

setMethod("[", c("GDset", "missing", "ANY"),
          function(x, i, j, ..., drop = FALSE){
            new("GDset", annot = x@annot, 
                dat = x@dat[ , j, drop=FALSE],
                pheno = x@pheno[j, , drop=FALSE], 
                platform = x@platform)  
          })

setMethod("[", c("GDset", "ANY", "missing"),
          function(x, i, j, ..., drop = FALSE){
            new("GDset", annot = x@annot[i], 
                dat = x@dat[i, , drop=FALSE],
                pheno = x@pheno[, , drop=FALSE], 
                platform = x@platform)  
          })

# Summaries ---------------------------------------------------------------

setMethod("show", "GDset", function(object) {
  cat("A GDset object \n")
  cat("Platform:", object@platform, "\n")
  cat("Data contains: \n")
  cat("  ", nrow(object@dat), "loci \n")
  cat("  ", ncol(object@dat), "samples \n")
  cat("With", ncol(object@pheno), "Covariates:\n")
  cat(colnames(object@pheno))
  
})

setMethod("dim", "GDset", function(x){
  c(loci = nrow(x@dat), samples = ncol(x@dat))
})



