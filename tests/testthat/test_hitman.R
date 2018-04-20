context("hitman")

test_that("hitman", {
  hm <- hitman(E=rnorm(length(pheno.v)), M=M, Y=pheno.v)
  expect_lt(mean(hm$comb.p < 0.05), 0.1)
  hm2 <- hitman(E=rnorm(length(pheno.v)), M=M, Y=pheno.v, covar = covar)
  expect_equal(mean(hm2$comb.p==hm$comb.p), 0)
})