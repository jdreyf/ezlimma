context("limma_dep")

test_that("bad Y", {
  expect_error(limma_dep(M, Y = rep("a", ncol(M))))
  expect_error(limma_dep(M, Y = rep(1, ncol(M))))
  expect_error(limma_dep(M, Y = c(1:(ncol(M)-1), NA)))
})

test_that("binary Y", {
  ld <- limma_dep(object=M, Y=y)
  expect_equal(ld$p, eztt[rownames(ld), "Last3vsFirst3.p"])
  ld2 <- limma_dep(object=M, Y=y, covar=covar)
  expect_equal(mean(ld2$p == eztt[rownames(ld2), "Last3vsFirst3.p"]), 0)
})

# lc <- limma_cor(M, phenotype = pheno.v) in helper
test_that("continuous Y", {
  ld <- limma_dep(object=M, Y = pheno.v)
  expect_equal(ld$p, lc$p)
  
  ld2 <- limma_dep(M, Y = pheno.v, covar=covar)
  expect_equal(mean(ld2$p == lc$p), 0)
})

test_that("matrix Y", {
  #compare to limmaF
  y3 <- cbind(y, pheno.v)
  rownames(y3) <- colnames(M)
  ld3 <- limma_dep(object=M, Y=y3)
  expect_equal(colnames(ld3)[1], "F")
  
  des <- model.matrix(~., data=data.frame(y3))
  lf <- limmaF(object = M, design = des, coef = 2:3)
  expect_equal(ld3, lf)
  #diff results without intercept
  LFnoInt <- limmaF(object = M, design = des[,-1])
  expect_equal(mean(ld3==LFnoInt), 0)
  
  #matrix Y with covar
  ld4 <- limma_dep(object=M, Y=y3, covar=covar)
  expect_equal(mean(ld3$p == ld4$p), 0)
})