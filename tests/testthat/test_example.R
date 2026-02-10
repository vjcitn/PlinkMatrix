
test_that("example is as expected", {
 ex = example_PlinkMatrix()
 expect_true(all(dim(ex) == c(445L, 367759L)))
})

#test_that("subset RangedPlink functions", {
# example(RangedPlink) # creates lit and myGR
# expect_true(max(IRanges::start((slot(IRanges::subsetByOverlaps(lit, myGR), "ranges")))) < 500000)
#})

test_that("subsetting of example with as_RSE is successful", {
 ex = example_PlinkMatrix(as_RSE=TRUE)
 GenomeInfoDb::seqlevelsStyle(rowRanges(ex)) = "Ensembl"
 simp = GenomicRanges::GRanges("18:1-200000")
 rs = subsetByOverlaps(ex, simp)
 expect_true(nrow(rs) == 835L)
})
