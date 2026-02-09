
test_that("example is as expected", {
 ex = example_PlinkMatrix()
 expect_true(all(dim(ex) == c(445L, 367759L)))
})

#test_that("subset RangedPlink functions", {
# example(RangedPlink) # creates lit and myGR
# expect_true(max(IRanges::start((slot(IRanges::subsetByOverlaps(lit, myGR), "ranges")))) < 500000)
#})


