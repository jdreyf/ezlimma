context("hitman")

#create associated phenotype, to avoid hitman warning about weak assoc
pheno.v <- setNames(rnorm(ncol(M)), nm=colnames(M))
pheno.v[1:3] <- pheno.v[1:3]-3
ee <- pheno.v + rnorm(length(pheno.v), sd=0.1)

test_that("E numeric", {
  expect_message(hm <- hitman(E=ee, M=M, Y=pheno.v, verbose = TRUE))
  expect_lt(mean(hm$EMY.p < 0.05), 0.1)
  
  covar.tmp <- rnorm(length(pheno.v))
  hm2 <- hitman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$EMY.p==hm2[rownames(hm), "EMY.p"]), 0.01)
  
  #no variance
  expect_error(hitman(E=numeric(length(pheno.v)), M=M, Y=pheno.v))
  
  ee2 <- rnorm(length(pheno.v), sd=0.1)
  expect_warning(hm3 <- hitman(E=ee2, M=M, Y=pheno.v))
  expect_lt(mean(hm3$EMY.p < 0.05), 0.1)
})

test_that("E binary", {
  expect_message(hm <- hitman(E=grp, M=M, Y=pheno.v, verbose = TRUE))
  expect_lt(mean(hm$EMY.p < 0.05), 0.2)
  
  #same as factor but not numeric
  hm2 <- hitman(E=factor(grp), M=M, Y=pheno.v)
  expect_equal(hm$EMY.p, hm2[rownames(hm), "EMY.p"])
  
  covar.tmp <- rnorm(length(pheno.v))
  expect_warning(hm3 <- hitman(E=grp, M=M, Y=pheno.v, covariates=covar.tmp))
  expect_lte(mean(hm$EMY.p==hm3[rownames(hm), "EMY.p"]), 0.01)
  
  y <- rep(1:3, times=2)
  limma_dep(object=y, Y=grp, prefix="EY")
  expect_warning(hm4 <- hitman(E=grp, M=M, Y=rep(1:3, times=2)))
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
  
  expect_equal(hm["gene1", "MY.p"], hm["gene1", "MY_dir.p"])
  expect_equal(hm["gene1", "EM.p"], hm["gene1", "EM_dir.p"])
  expect_equal(hm["gene1", "EMY.p"], max(hm["gene1", "EM.p"], hm["gene1", "MY.p"])^2)
})

test_that("NAs", {
  expect_error(hitman(E=grp, M=M, Y=pheno2))
  
  grp2 <- grp
  grp2[1] <- NA
  expect_error(hitman(E=grp2, M=M, Y=pheno.v))
  
  M[1, 1] <- NA
  expect_silent(hitman(E=ee, M=M, Y=pheno.v))
})