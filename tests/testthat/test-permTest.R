context("Testing permTest function")

# create components to test GDset creation
annot <- GRanges(seqnames = Rle(rep(c('chr1', 'chr8'), each=5)),
                 ranges = IRanges(start = 20*(1:10), end = 20*(1:10)+10),
                 entrez.id = rep(c('HOX', 'FOX'), each=5))

names(annot) <- paste0('probe', 1:10)

expData <- data.frame(sample1 = rnorm(10) + 1:10, sample2 = rnorm(10) + 1:10,
                      sample3 = rnorm(10) + 1:10, sample4 = rnorm(10) + 1:10,
                      sample5 = rnorm(10), sample6 = rnorm(10),
                      sample7 = rnorm(10), sample8 = rnorm(10))

row.names(expData) <- names(annot)

pData <- data.frame(batch = rep(1:4, 2), 
                    region = rep(c("DFC", "CBC"), each=4))

rownames(pData) <- colnames(expData)

GDset1 <- GDset(annot = annot, dat = expData, 
                pheno = pData, platform = "microarray")

annot2 <- annot[1:7]
expData2 <- expData[1:7, ]
expData2 <- expData2 + rnorm(length(expData2))*3
GDset2 <- GDset(annot = annot2, dat = expData2,
                pheno = pData, platform = 'methy')

GDIset.test <- GDIset(GDset1, GDset2)


permTest(GDIset.test, n.perm = 1e4)
# 
# test_that("'permTest' returns numeric", {
#   expect_that(class(permTest(GDIset.test)), equals('numeric'))
# })



