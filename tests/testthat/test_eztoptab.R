context('eztoptab')


fit <- eBayes(lmFit(M))
test <- ezlimma:::eztoptab(fit,adjust.method = 'none', prefix = 'abc')


test_that("column names are adjusted to adjust.method and prefixes added", {
  
  
  expect_equal(FALSE, any(grepl('adj\\.P\\.Val', colnames(test))))
  expect_equal(TRUE,any(grepl('abc',colnames(test))))
})
