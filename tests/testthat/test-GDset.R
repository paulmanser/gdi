context("Testing methods for GDset")

# create components to test GDset creation
annot <- GRanges(seqnames = Rle(c('chr1', 'chr8')),
                 ranges = IRanges(start = 1:2, end = 4:5),
                 entrez.id = c('HOX', 'FOX'))
names(annot) <- c('probe1', 'probe2')

expData <- data.frame(sample1 = 1:2, sample2 = 4:5,
                   sample3 = 8:9, sample4 = 0:1)
row.names(expData) <- c('probe1', 'probe2')

pData <- data.frame(batch = 1:4, region = rep(c("DFC", "CBC"), each=2))
rownames(pData) <- colnames(expData)

options(fftempdir = 'C:/Users/pmanser/Documents/gdi')
expData <- as.ffdf(expData)

GDset.test <- GDset(annot = annot, dat = expData, 
                    pheno = pData, platform = "microarray")

# GDset creation
test_that("GDset creation", {
  expect_that(class(GDset.test)[1], equals("GDset"))
})

test_that("Throw error on GDset mismatch", {
  expect_that(GDset(annot = annot, dat = expData[1, ], 
                    pheno = pData, platform = "microarray"), throws_error())
})

test_that("Throw error on GDset mismatch", {
  expect_that(GDset(annot = annot, dat = expData, 
                    pheno = pData[1, ], platform = "microarray"), throws_error())
})

test_that("Throw error on GDset mismatch", {
  expect_that(GDset(annot = annot[1], dat = expData, 
                    pheno = pData, platform = "microarray"), throws_error())
})

# accessors
test_that("'getAnnot' accessor", {
  expect_that(getAnnot(GDset.test), equals(annot))
})

test_that("'getDat' accessor", {
  expect_that(identical(getDat(GDset.test), expData), equals(TRUE))
})

test_that("'getPheno' accessor", {
  expect_that(getPheno(GDset.test), equals(pData))
})

test_that("'getPlatform' accessor", {
  expect_that(getPlatform(GDset.test), equals("microarray"))  
})

# subsetting
dim.test <- c(loci = 2, samples = 4)
test_that("Test dim()", {
  expect_that(dim(GDset.test), equals(dim.test))
})

dim.test <- c(loci = 1, samples = 4)
test_that("Subsetting GDset rows", {
  expect_that(dim(GDset.test[1, ]), equals(dim.test))
})

dim.test <- c(loci = 2, samples = 2)
test_that("Subsetting GDset cols", {
  expect_that(dim(GDset.test[ , 2:3]), equals(dim.test))
})

dim.test <- c(loci = 1, samples = 1)
test_that("Subsetting rows & cols", {
  expect_that(dim(GDset.test[1, 1]), equals(dim.test))
})

test_that("subsetting rows & cols by name", {
  expect_that(dim(GDset.test['probe1', 'sample1']), equals(dim.test))
})

