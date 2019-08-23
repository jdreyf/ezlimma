context('add z cols')

# if fit holds > 1 comparison, eztoptab yields F-stat
test_that("one comparison w/ prefix", {
  eztt1 <- ezlimma:::eztoptab(fit, prefix = 'abc', cols = c("t", "logFC", "P.Value"))
  eztt1.z <- add_zcols(eztt1)
  expect_true("abc.z" %in% colnames(eztt1.z))
  expect_lte(which(colnames(eztt1.z) == "abc.z"), which(colnames(eztt1.z) == "abc.t"))
  expect_true(all(abs(eztt1.z$abc.z) <= abs(eztt1.z$abc.t)))
})

test_that("one comparison w/o prefix", {
  eztt2 <- ezlimma:::eztoptab(fit, cols = c("t", "P.Value"))
  eztt2.z <- add_zcols(eztt2)
  expect_true("z" %in% colnames(eztt2.z))
  expect_lte(which(colnames(eztt2.z) == "z"), which(colnames(eztt2.z) == "t"))
  expect_true(all(abs(eztt2.z$z) <= abs(eztt2.z$t)))
})

test_that(">1 comparisons w/ prefix  for logFC", {
  eztt.z <- add_zcols(eztt, stat.suffix = "logFC")
  expect_true(all(c("First3.z", "Last3.z", "Last3vsFirst3.z") %in% colnames(eztt.z)))
  expect_true(all(sign(eztt.z$First3.z) <= sign(eztt.z$First3.logFC)))
})