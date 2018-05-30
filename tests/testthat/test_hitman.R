context("hitman")

test_that("hitman", {
  ee <- rnorm(length(pheno.v))
  hm <- hitman(E=ee, M=M, Y=pheno.v)
  expect_lt(mean(hm$comb.p < 0.05), 0.1)
  
  covar.tmp <- rnorm(length(pheno.v))
  hm2 <- hitman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$comb.p==hm2[rownames(hm), "comb.p"]), 0.01)
})