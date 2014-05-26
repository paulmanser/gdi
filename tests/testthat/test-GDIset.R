context("Testing methods for GDIset")

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


GDset1 <- GDset(annot = annot, dat = expData, 
                pheno = pData, platform = "microarray")

annot2 <- GRanges(seqnames = Rle(c('chr1', 'chr8', 'chr2')),
                  ranges = IRanges(start = 1:3, end = 4:6),
                  entrez.id = c('HOX', 'FOX', 'GABA'))

names(annot2) <- paste0('probe', 1:3)

expData2 <- data.frame(sample1 = 1:3, sample2 = 4:6,
                       sample3 = 8:10, sample4 = 0:2)
rownames(expData2) <- names(annot2)

GDset2 <- GDset(annot = annot2, dat = expData2,
                pheno = pData, platform = 'methy')

GDIset.test <- GDIset(GDset1, GDset2)

# creation
test_that("GDIset creation", {
  expect_that(class(GDIset.test)[1], equals('GDIset'))
})

test_that("GDIset mismatched samples", {
  expect_that(GDIset(GDset1[, 1:2], GDset2), throws_error())
})

#accessors
annot.test <- list(microarray = annot, methy = annot2)
test_that("'getAnnot' accessor", {
  expect_that(getAnnot(GDIset.test), equals(annot.test))
})

dat.test <- list(microarray = expData, methy = expData2)
test_that("'getDat' accessor", {
  expect_that(getDat(GDIset.test), equals(dat.test))
})

test_that("'getPheno' accessor", {
  expect_that(getPheno(GDIset.test), equals(pData))
})

plat.test <- list(set1 = "microarray", set2 = "methy")
test_that("'getPlatform' accessor", {
  expect_that(getPlatform(GDIset.test), equals(plat.test))  
})

test_that("getSet accessor for first slot", {
  expect_that(getSet(GDIset.test, 1), equals(GDset1))
})

test_that("getSet accessor for second slot", {
  expect_that(getSet(GDIset.test, 2), equals(GDset2))
})

# subsetting
hox.dim <- list(microarray=c(loci=1, samples=4),
                methy=c(loci=1, samples=4))
test_that("Subset GDIset by gene symbol", {
  expect_that(dim(GDIset.test['HOX', ]), equals(hox.dim))
})

hox.dim <- list(microarray=c(loci=2, samples=2),
                methy=c(loci=3, samples=2))
test_that("Subset GDIset by gene symbol", {
  expect_that(dim(GDIset.test[ , 2:3]), equals(hox.dim))
})

hox.dim <- list(microarray=c(loci=1, samples=1),
                methy=c(loci=1, samples=1))
test_that("Subset GDIset by gene symbol", {
  expect_that(dim(GDIset.test['HOX', 1]), equals(hox.dim))
})

hox.dim <- list(microarray=c(loci=1, samples=4),
                methy=c(loci=1, samples=4))
test_that("Subset GDIset by gene symbol", {
  expect_that(GDIset.test[5, ], throws_error())
})

# consolidate
test_that("Consolidate GDIset returns GDIset", {
  expect_that(class(consolidate(GDIset.test))[1], equals("GDIset"))
})

x <- consolidate(GDIset.test)
test_that("Consolidate returns correctly", {
  expect_that(identical(unique(x@set1@annot$entrez.id),
                        unique(x@set2@annot$entrez.id)),
    is_true())
})



