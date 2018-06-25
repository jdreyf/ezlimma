context('eztoptab')

fit <- limma::eBayes(lmFit(M))

test_that("adjust.method none", {
  eztt.none <- ezlimma:::eztoptab(fit, adjust.method = 'none', prefix = 'abc')
  expect_false(any(grepl('adj\\.P\\.Val', colnames(eztt.none))))
  expect_true(all(grepl('abc', colnames(eztt.none))))
  expect_equal(eztt.none$abc.p, eztt.none$abc.none)
})

test_that("adjust.method holm", {
  eztt.holm <- ezlimma:::eztoptab(fit, adjust.method = 'holm', prefix = 'a-b')
  expect_false(any(grepl('adj\\.P\\.Val', colnames(eztt.holm))))
  expect_true(all(grepl('a-b', colnames(eztt.holm), fixed = TRUE)))
  expect_true(all(eztt.holm$`a-b.p` <= eztt.holm$`a-b.holm`))
})
