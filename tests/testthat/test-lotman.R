context("lotman")

# create associated phenotype, to avoid lotman warning about weak assoc
set.seed(0)
pheno.v <- setNames(rnorm(ncol(M)), nm=colnames(M))
pheno.v[1:3] <- pheno.v[1:3]-3
ee <- pheno.v + rnorm(length(pheno.v), sd=0.1)
grp2 <- batch2design(grp)[,1]
names(grp2) <- colnames(M)
covar.tmp <- rnorm(length(pheno.v))
# E=ee; Y=pheno.v; covariates=covar.tmp
# E=grp.tmp;M=M;Y=pheno.v;covariates=covar.tmp

test_that("E numeric", {
  hm <- lotman(E=ee, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.1)
  
  hm2 <- lotman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$EMY.p==hm2[rownames(hm), "EMY.p"]), 0.01)
  
  #no variance
  expect_error(lotman(E=numeric(length(pheno.v)), M=M, Y=pheno.v))
  
  ee2 <- rnorm(length(pheno.v), sd=0.1)
  expect_warning(hm3 <- lotman(E=ee2, M=M, Y=pheno.v))
  expect_lt(mean(hm3$EMY.p < 0.05), 0.1)
})

test_that("E binary", {
  hm <- lotman(E=grp2, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.2)
  
  #EMY.p indep of parametrization
  grp3 <- as.numeric(grp==grp[1])
  hm2 <- lotman(E=grp3, M=M, Y=pheno.v)
  expect_equal(hm$EMY.p, hm2$EMY.p)

  covar2 <- rnorm(length(pheno.v))
  expect_warning(hm3 <- lotman(E=grp2, M=M[1,], Y=pheno.v, covariates=covar2)) # no assoc warning
  expect_true(hm["gene1", "EMY.p"] != hm3[1, "EMY.p"])
  
  y <- rep(1:3, times=2)
  expect_warning(hm4 <- lotman(E=grp2, M=M, Y=rep(1:3, times=2)))
})

test_that("E nominal --> design", {
  grp.tmp <- batch2design(rep(letters[1:2], each=3))[,1]
  names(grp.tmp) <- colnames(M)
  
  hm <- lotman(E=grp.tmp, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.2)
  
  # warning: essentially perfect fit: summary may be unreliable
  expect_warning(hm3 <- lotman(E=grp.tmp, M=M, Y=pheno.v, covariates=covar.tmp))
  expect_lte(mean(hm$EMY.p == hm3[rownames(hm), "EMY.p"]), 0.01)
  
  expect_error(lotman(E=rep("a", length(pheno.v)), M=M, Y=pheno.v))
  expect_error(lotman(E=c(rep("a", length(pheno.v)-1), NA), M=M, Y=pheno.v))
  expect_error(lotman(E=rep(NA, length(pheno.v)-1), M=M, Y=pheno.v))
})

test_that("M df", {
  mdf <- data.frame(M)
  hm <- lotman(E=ee, M=mdf, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.1)
})

test_that("gene1", {
  hm <- lotman(E=grp2, M=M, Y=M[1,])
  expect_equal(rownames(hm)[1], "gene1")
  
  expect_equal(hm["gene1", "MY.p"], hm["gene1", "MY_dir.p"])
  expect_equal(hm["gene1", "EM.p"], hm["gene1", "EM_dir.p"])
  expect_equal(hm["gene1", "EMY.p"], max(hm["gene1", "EM.p"], hm["gene1", "MY.p"])^2)
})

test_that("NAs", {
  expect_error(lotman(E=grp2, M=M, Y=pheno2))
  
  grp2[1] <- NA
  expect_error(lotman(E=grp2, M=M, Y=pheno.v))
  
  M[1, 1] <- NA
  expect_silent(lotman(E=ee, M=M, Y=pheno.v))
})

test_that("testing a gene is independent of other genes", {
  hm1 <- lotman(E=ee, M=M[1,], Y=pheno.v, covariates=covar.tmp)
  hm2 <- lotman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_equal(hm1[1, "EMY.p"], hm2["gene1", "EMY.p"])
})

test_that("consistent & inconsistent", {
  n <- 10
  sigma <- 0.25
  E <- rep(0:1, each=n)
  ey <- rnorm(n=2*n, sd=sigma)
  Y <- E+ey
  eps <- rnorm(n=2*n, sd=sigma)
  em.ic <- -ey+eps
  em.c <- ey+eps
  M <- rbind(ics=E+em.ic, cs=E+em.c)
  lm <- lotman(E=E, M=M, Y=Y)
  expect_lt(lm["cs", "EMY.p"], 0.01)
  expect_gt(lm["ics", "EMY.p"], 0.9)
})

# takes a few sec -- worth it.
test_that("barfield", {
  prop.sig.mat <- ezlimma:::sim_barfield(med.fcn = lotman, b1t2.v=c(0, 0.39), nsim = 50)
  expect_lte(prop.sig.mat[1, 1], 0.03)
  expect_gte(prop.sig.mat[2, 2], 0.6)
})