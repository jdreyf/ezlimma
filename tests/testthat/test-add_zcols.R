context('add z cols')

# if fit holds > 1 comparison, eztoptab yields F-stat
test_that("one comparison", {
  eztt1 <- ezlimma:::eztoptab(fit, prefix = 'abc', cols = c("t", "logFC", "P.Value"))
  eztt2 <- add_zcols(eztt1, fit)
  expect_true("abc.z" %in% colnames(eztt2))
  expect_lte(which(colnames(eztt2) == "abc.z"), which(colnames(eztt2) == "abc.t"))
  expect_true(all(abs(eztt2$abc.z) <= abs(eztt2$abc.t)))
})