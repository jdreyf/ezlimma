context("hitman")

test_that("E numeric", {
  ee <- rnorm(length(pheno.v))
  expect_message(hm <- hitman(E=ee, M=M, Y=pheno.v, verbose = TRUE))
  expect_lt(mean(hm$EMY.p < 0.05), 0.1)
  
  covar.tmp <- rnorm(length(pheno.v))
  hm2 <- hitman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$EMY.p==hm2[rownames(hm), "EMY.p"]), 0.01)
  
  expect_error(hitman(E=numeric(length(pheno.v)), M=M, Y=pheno.v))
})

test_that("E binary", {
  expect_message(hm <- hitman(E=grp, M=M, Y=pheno.v, verbose = TRUE))
  expect_lt(mean(hm$EMY.p < 0.05), 0.2)
  
  #same as factor but not numeric
  hm2 <- hitman(E=factor(grp), M=M, Y=pheno.v)
  expect_equal(hm$EMY.p, hm2[rownames(hm), "EMY.p"])
  
  covar.tmp <- rnorm(length(pheno.v))
  hm3 <- hitman(E=grp, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$EMY.p==hm3[rownames(hm), "EMY.p"]), 0.01)
  
  #try to get ey.t=0, but EY.t~=1e-16
  # y <- rep(1:3, times=2)
  # limma_dep(object=y, Y=grp, prefix="EY")
  # hm4 <- hitman(E=grp, M=M, Y=rep(1:3, times=2))
})

test_that("E nominal", {
  grp.tmp <- rep(letters[1:3], each=2)
  
  hm <- hitman(E=grp.tmp, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.2)
  
  #same as factor but not numeric
  hm2 <- hitman(E=factor(grp.tmp), M=M, Y=pheno.v)
  expect_equal(hm$EMY.p, hm2[rownames(hm), "EMY.p"])
  
  covar.tmp <- rnorm(length(pheno.v))
  hm3 <- hitman(E=grp.tmp, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$EMY.p==hm3[rownames(hm), "EMY.p"]), 0.01)
  
  expect_error(hitman(E=rep("a", length(pheno.v)), M=M, Y=pheno.v))
  expect_error(hitman(E=c(rep("a", length(pheno.v)-1), NA), M=M, Y=pheno.v))
  expect_error(hitman(E=rep(NA, length(pheno.v)-1), M=M, Y=pheno.v))
})

test_that("gene1", {
  hm <- hitman(E=grp, M=M, Y=M[1,])
  expect_equal(rownames(hm)[1], "gene1")
  
  expect_equal(hm["gene1", "MY.p"], hm["gene1", "MY2.p"])
  expect_equal(hm["gene1", "EM.p"], hm["gene1", "EM2.p"])
  expect_equal(hm["gene1", "EMY.p"], max(hm["gene1", "EM.p"], hm["gene1", "MY.p"])^2)
})