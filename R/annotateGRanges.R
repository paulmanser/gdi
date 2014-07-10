

annotateGRanges <- function(gr){
  
  if (!is(gr, 'GenomicRanges'))
    stop("'gr' must be a 'GenomicRanges' object")
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  # get exon to gene_id key
  exons <- unlist(exonsBy(txdb, "gene"))
  exon2gene <- names(exons)
  names(exon2gene) <- exons$exon_id
  
  # get tx to gene_id key
  tx <- unlist(transcriptsBy(txdb, "gene"))
  tx2gene <- names(tx)
  names(tx2gene) <- tx$tx_id
  
  # 3' UTRs ---------------------------------------------
  utr3 <- unlist(threeUTRsByTranscript(txdb))
  gene.ids <- exon2gene[as.character(utr3$exon_id)]
  utr3$gene_id <- gene.ids
  utr3 <- utr3[-which(is.na(gene.ids))]
  gene.ids <- na.omit(gene.ids)
  utr.by.gene <- unlist(reduce(split(utr3, as.character(gene.ids))))
  utr.by.gene$gene_id <- names(utr.by.gene)
  utr.by.gene$location <- "3UTR"
  utr3 <- utr.by.gene
  
  # 5' UTRs ---------------------------------------------
  utr5 <- unlist(fiveUTRsByTranscript(txdb))
  gene.ids <- exon2gene[as.character(utr5$exon_id)]
  utr5$gene_id <- gene.ids
  utr5 <- utr5[-which(is.na(gene.ids))]
  gene.ids <- na.omit(gene.ids)
  utr.by.gene <- unlist(reduce(split(utr5, as.character(gene.ids))))
  utr.by.gene$gene_id <- names(utr.by.gene)
  utr.by.gene$location <- "5UTR"
  utr5 <- utr.by.gene
   
  # promoters -------------------------------------------
  proms <- promoters(txdb, upstream=3000, downstream=100)
  gene.ids <- tx2gene[as.character(proms$tx_id)]
  proms$gene_id <- gene.ids
  proms <- proms[-which(is.na(gene.ids))]
  gene.ids <- na.omit(gene.ids)
  prom.by.gene <- unlist(reduce(split(proms, as.character(gene.ids))))
  prom.by.gene$gene_id <- names(prom.by.gene)
  prom.by.gene$location <- "Promoter"
  proms <- prom.by.gene
   
  # coding ----------------------------------------------
  cds <- unlist(reduce(cdsBy(txdb, "gene")))
  cds$gene_id <- names(cds)
  cds$location <- "Coding"
   
  # intron ----------------------------------------------
  introns <- unlist(intronsByTranscript(txdb))
  gene.ids <- tx2gene[as.character(names(introns))]
  introns$gene_id <- gene.ids
  introns <- introns[-which(is.na(gene.ids))]
  gene.ids <- na.omit(gene.ids)
  intron.by.gene <- unlist(reduce(split(introns, as.character(gene.ids))))
  intron.by.gene$gene_id <- names(intron.by.gene)
  intron.by.gene$location <- "Intron"
  introns <- intron.by.gene
  
  hg19.annot.gr <- c(utr3, utr5, proms, cds, introns)
  
  
  # overlap with annotation ------------------------------------
  
  annot.overlaps <- findOverlaps(hg19.annot.gr, gr)
  
  q.split <- split(queryHits(annot.overlaps),
                   subjectHits(annot.overlaps))
  
  annotatedRanges <- gr
  annotatedRanges$location <- annotatedRanges$gene_id <- NA
  
  chr.gene.ids <- tapply(hg19.annot.gr$gene_id[queryHits(annot.overlaps)],
                         factor(as.character(subjectHits(annot.overlaps))),
                         function(x) paste0(x, collapse=';'))
  
  chr.locs <- tapply(hg19.annot.gr$location[queryHits(annot.overlaps)],
                     factor(as.character(subjectHits(annot.overlaps))),
                     function(x) paste0(x, collapse=';'))
  
  chr.strand <- tapply(as.character(strand(hg19.annot.gr))[queryHits(annot.overlaps)],
                       factor(as.character(subjectHits(annot.overlaps))),
                       function(x) x[1])
  
  annotatedRanges$gene_id[as.numeric(names(chr.gene.ids))] <- chr.gene.ids
  annotatedRanges$location[as.numeric(names(chr.locs))] <- chr.locs
  strand(annotatedRanges)[as.numeric(names(chr.strand))] <- chr.strand
  
  annotatedRanges
}

