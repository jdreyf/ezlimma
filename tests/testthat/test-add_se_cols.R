context('add se cols')

# if fit holds > 1 comparison, eztoptab yields F-stat
test_that("one comparison w/ prefix", {
  eztt1 <- ezlimma:::eztoptab(fit, prefix = 'abc', cols = c("t", "logFC", "P.Value"))
  eztt1.SE <- add_se_cols(eztt1)
  expect_true("abc.SE" %in% colnames(eztt1.SE))
  expect_gt(which(colnames(eztt1.SE) == "abc.SE"), which(colnames(eztt1.SE) == "abc.t"))
  expect_true(all(eztt1.SE$abc.SE == eztt1.SE$abc.logFC/eztt1.SE$abc.t))
})

test_that("one comparison w/o prefix", {
  eztt2 <- ezlimma:::eztoptab(fit, cols = c("t", "logFC", "P.Value"))
  eztt2.SE <- add_se_cols(eztt2)
  expect_true("SE" %in% colnames(eztt2.SE))
  expect_gt(which(colnames(eztt2.SE) == "SE"), which(colnames(eztt2.SE) == "t"))
  expect_true(all(eztt2.SE$SE == eztt2.SE$logFC/eztt2.SE$t))
})

test_that(">1 comparisons w/ prefix  for logFC", {
  eztt.t <- limma_contrasts(object=M, grp = grp, contrast.v = contr.v, cols = c("t", "logFC", "P.Value"))
  eztt.SE <- add_se_cols(eztt.t)
  expect_true(all(c("First3.SE", "Last3.SE", "Last3vsFirst3.SE") %in% colnames(eztt.SE)))
  expect_true(all(eztt.t$SE == eztt.t$logFC/eztt.t$t))
})