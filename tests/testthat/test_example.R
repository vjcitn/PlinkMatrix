
test_that("example is as expected", {
 ex = example_PlinkArray()
 expect_true(all(dim(ex) == c(445L, 367759L)))
})
