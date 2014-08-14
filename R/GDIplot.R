

GDIplot <- function(object, covariate=NULL, covariate.type=NULL){
  
  if (!is(object, 'GDIset')) stop("object must be a 'GDIset'")
  
  if (is.null(covariate) | !any(covariate %in% colnames(getPheno(object))))
    stop(paste0("Specify a 'covariate' name from GDI object 'pheno' slot.
                Available covariates are: ", paste(colnames(getPheno(object)), collapse=', ')))
  
  
  if (is.null(covariate.type) | !any(covariate.type %in% c("continuous", "categorical")))
    stop("Specify a covariate type. 'covariate.type' must be either 'continuous' or 'categorical'")
  
  # create layout ----------------------------------------------
  layout.mat <- matrix(2, 5, 5)
  layout.mat[, 1] <- 1
  layout(layout.mat)
  
  # create gene model ------------------------------------------
  if (as.character(strand(object@set1@annot)[1]) == '-'){
    strand.top <- "5'"
    strand.bot <- "3'"
  } else {
    strand.top <- "3'"
    strand.bot <- "5'"
  }
  
  par(mar = c(4, 3, 4, 0))
  plot(0, 0, type = 'n', axes = FALSE, xlab=strand.bot, ylab='Distance (kb)',
       ylim = ylims, xlim = c(0, 1), main = strand.top)
  
  axis(2, at = seq(ylims[1], ylims[2], length.out = 10),
       labels = round(seq(ylims[1], ylims[2], length.out = 10)/1e3, 1))
  
  # add gene model on left -----------------------------------
  # map gene symbol to entrez gene
  
  exon.gr <- unlist(reduce(exons[eg.ids]))
  cds.gr <- unlist(reduce(cds[eg.ids]))
    
  gene.range <- range(exon.gr)
  gene.range <- c(start(gene.range), end(gene.range))
  segments(.5, gene.range[1], .5, gene.range[2])
  
  # add exon which include UTR
  for(jj in 1:length(exon.gr)){
    gr <- c(start(exon.gr[jj]), end(exon.gr[jj]))
    rect(.45, min(gr), .55, max(gr), col = 1)
  }
  
  # add cds
  for(jj in 1:length(cds.gr)){
    gr <- c(start(cds.gr[jj]), end(cds.gr[jj]))
    rect(.4, min(gr), .6, max(gr), col = 1)
  }
  
  
  
  
  
}