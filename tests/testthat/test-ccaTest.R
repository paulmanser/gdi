context("Testing ccaTest function")
options("fftempdir"=getwd())
# create components to test GDset creation
annot <- GRanges(seqnames = Rle(rep(c('chr1', 'chr8'), each=5)),
                 ranges = IRanges(start = 1:10, end = 4*(1:10)+10),
                 entrez.id = rep(c('HOX', 'FOX'), each=5))

names(annot) <- paste0('probe', 1:10)

expData <- data.frame(sample1 = rnorm(10) + 1:10, sample2 = rnorm(10) + 1:10,
                      sample3 = rnorm(10) + 1:10, sample4 = rnorm(10) + 1:10,
                      sample5 = rnorm(10), sample6 = rnorm(10),
                      sample7 = rnorm(10), sample8 = rnorm(10))

row.names(expData) <- names(annot)
expData <- as.ffdf(expData)

pData <- data.frame(batch = rep(1:4, 2), 
                    region = rep(c("DFC", "CBC"), each=4))
pData <- pData

rownames(pData) <- colnames(expData)

GDset1 <- GDset(annot = annot, dat = expData, 
                pheno = pData, platform = "microarray")

annot2 <- annot[1:7]
expData2 <- expData[1:7, ]
expData2 <- expData2 + rnorm(length(expData2))*3
expData2 <- as.ffdf(expData2)

GDset2 <- GDset(annot = annot2, dat = expData2,
                pheno = pData, platform = 'methy')

GDIset.test <- GDIset(GDset1, GDset2)


ccaTest(GDIset.test)




test_that("'ccaTest' returns list", {
  expect_that(class(ccaTest(GDIset.test)), equals("list"))
})







