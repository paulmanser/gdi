

annotateGRanges <- function(gr, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene){
  
  if (!is(gr, 'GenomicRanges'))
    stop("'gr' must be a 'GenomicRanges' object")
  
  mapped.loci <- locateVariants(gr, txdb, AllVariants())
  unique.overlaps <- findOverlaps(unique(mapped.loci), mapped.loci)
  
  annot <- matrix(NA, nr = length(gr), nc = 2)
  colnames(annot) <- c('location', 'entrez.id')
  sh <- subjectHits(unique.overlaps)
  qh <- queryHits(unique.overlaps)
  
  for(ii in 1:length(gr)){
    location <- paste0(mapped.loci$LOCATION[sh[which(qh == ii)]], collapse=';')
    entrez.id <- paste0(mapped.loci$GENEID[sh[which(qh == ii)]], collapse=';')
    annot[ii, ] <- c(location, entrez.id)
  }
  
  gr$location <- annot[, 1]
  gr$entrez.id <- annot[, 2]
  gr
}

