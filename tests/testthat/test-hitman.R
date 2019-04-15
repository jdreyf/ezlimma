context("hitman")

# create associated phenotype, to avoid hitman warning about weak assoc
set.seed(0)
pheno.v <- setNames(rnorm(ncol(M)), nm=colnames(M))
pheno.v[1:3] <- pheno.v[1:3]-3
ee <- pheno.v + rnorm(length(pheno.v), sd=0.1)
grp2 <- batch2design(grp)[,1]
names(grp2) <- colnames(M)
covar.tmp <- rnorm(length(pheno.v))

test_that("E numeric", {
  hm <- hitman(E=ee, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.1)
  
  hm2 <- hitman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$EMY.p==hm2[rownames(hm), "EMY.p"]), 0.01)
  
  #no variance
  expect_error(hitman(E=numeric(length(pheno.v)), M=M, Y=pheno.v))
  
  ee2 <- rnorm(length(pheno.v), sd=0.1)
  expect_warning(hm3 <- hitman(E=ee2, M=M, Y=pheno.v))
  expect_lt(mean(hm3$EMY.p < 0.05), 0.1)
})

test_that("E binary", {
  hm <- hitman(E=grp2, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.2)
  
  #EMY.p indep of parametrization
  grp3 <- as.numeric(grp==grp[1])
  hm2 <- hitman(E=grp3, M=M, Y=pheno.v)
  expect_equal(hm$EMY.p, hm2$EMY.p)
  
  expect_warning(hm3 <- hitman(E=grp2, M=M, Y=pheno.v, covariates=covar.tmp))
  expect_lte(mean(hm$EMY.p == hm3[rownames(hm), "EMY.p"]), 0.01)
  
  y <- rep(1:3, times=2)
  expect_warning(hm4 <- hitman(E=grp2, M=M, Y=rep(1:3, times=2)))
})

test_that("E nominal --> design", {
  grp.tmp <- batch2design(rep(letters[1:2], each=3))[,1]
  names(grp.tmp) <- colnames(M)
  
  hm <- hitman(E=grp.tmp, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.2)
  
  set.seed(0)
  covar.tmp <- rnorm(length(pheno.v))
  # warning: essentially perfect fit: summary may be unreliable
  expect_warning(hm3 <- hitman(E=grp.tmp, M=M, Y=pheno.v, covariates=covar.tmp))
  expect_lte(mean(hm$EMY.p == hm3[rownames(hm), "EMY.p"]), 0.01)
  
  expect_error(hitman(E=rep("a", length(pheno.v)), M=M, Y=pheno.v))
  expect_error(hitman(E=c(rep("a", length(pheno.v)-1), NA), M=M, Y=pheno.v))
  expect_error(hitman(E=rep(NA, length(pheno.v)-1), M=M, Y=pheno.v))
})

test_that("M df", {
  mdf <- data.frame(M)
  hm <- hitman(E=ee, M=mdf, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.1)
})

test_that("gene1", {
  hm <- hitman(E=grp2, M=M, Y=M[1,])
  expect_equal(rownames(hm)[1], "gene1")
  
  expect_equal(hm["gene1", "MY.p"], hm["gene1", "MY_dir.p"])
  expect_equal(hm["gene1", "EM.p"], hm["gene1", "EM_dir.p"])
  expect_equal(hm["gene1", "EMY.p"], max(hm["gene1", "EM.p"], hm["gene1", "MY.p"])^2)
})

test_that("NAs", {
  expect_error(hitman(E=grp2, M=M, Y=pheno2))
  
  grp2[1] <- NA
  expect_error(hitman(E=grp2, M=M, Y=pheno.v))
  
  M[1, 1] <- NA
  expect_silent(hitman(E=ee, M=M, Y=pheno.v))
})

# takes a few sec -- worth it.
test_that("barfield", {
  prop.sig.mat <- ezlimma:::sim_barfield(med.fcn = hitman, b1t2.v=c(0, 0.39), nsim = 50, ngene = 9)
  expect_lte(prop.sig.mat[1, 1], 0.03)
  expect_gte(prop.sig.mat[2, 2], 0.6)
})