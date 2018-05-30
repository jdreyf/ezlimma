context("hitman")

test_that("E numeric", {
  ee <- rnorm(length(pheno.v))
  hm <- hitman(E=ee, M=M, Y=pheno.v)
  expect_lt(mean(hm$comb.p < 0.05), 0.1)
  
  covar.tmp <- rnorm(length(pheno.v))
  hm2 <- hitman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$comb.p==hm2[rownames(hm), "comb.p"]), 0.01)
})

test_that("E binary", {
  hm <- hitman(E=grp, M=M, Y=pheno.v)
  expect_lt(mean(hm$comb.p < 0.05), 0.2)
  
  #same as factor but not numeric
  hm2 <- hitman(E=factor(grp), M=M, Y=pheno.v)
  expect_equal(hm$comb.p, hm2[rownames(hm), "comb.p"])
  
  covar.tmp <- rnorm(length(pheno.v))
  hm3 <- hitman(E=grp, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$comb.p==hm3[rownames(hm), "comb.p"]), 0.01)
})

test_that("E nominal", {
  grp.tmp <- rep(letters[1:3], each=2)
  
  hm <- hitman(E=grp.tmp, M=M, Y=pheno.v)
  expect_lt(mean(hm$comb.p < 0.05), 0.2)
  
  #same as factor but not numeric
  hm2 <- hitman(E=factor(grp.tmp), M=M, Y=pheno.v)
  expect_equal(hm$comb.p, hm2[rownames(hm), "comb.p"])
  
  covar.tmp <- rnorm(length(pheno.v))
  hm3 <- hitman(E=grp.tmp, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$comb.p==hm3[rownames(hm), "comb.p"]), 0.01)
})

test_that("gene1", {
  hm <- hitman(E=grp, M=M, Y=M[1,])
  expect_equal(rownames(hm)[1], "gene1")
  
  expect_equal(hm["gene1", "MY.p"], hm["gene1", "MY2.p"])
  expect_equal(hm["gene1", "EM.p"], hm["gene1", "EM2.p"])
  expect_equal(hm["gene1", "comb.p"], max(hm["gene1", "EM.p"], hm["gene1", "MY.p"])^2)
})