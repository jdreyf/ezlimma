context("limma_dep")

test_that("Y", {
  expect_error(limma_dep(M, Y = rep("a", ncol(M))))
  expect_error(limma_dep(M, Y = rep(1, ncol(M))))
  expect_error(limma_dep(M, Y = c(1:(ncol(M)-1), NA)))
  
  ld <- limma_dep(object=M, Y=y)
  expect_equal(ld$p, eztt[rownames(ld), "Last3vsFirst3.p"])
  ld2 <- limma_dep(object=M, Y=y, covar=covar)
  expect_equal(mean(ld2$p == eztt[rownames(ld2), "Last3vsFirst3.p"]), 0)
})

# lc <- limma_cor(M, phenotype = pheno.v) in helper
test_that("covariate", {
  ld <- limma_dep(object=M, Y = pheno.v)
  expect_equal(ld$p, lc$p)
  
  ld2 <- limma_dep(M, Y = pheno.v, covar=covar)
  expect_equal(mean(ld2$p == lc$p), 0)
  
  #matrix Y
  y3 <- cbind(y, pheno.v)
  rownames(y3) <- colnames(M)
  ld3 <- limma_dep(object=M, Y=y3)
  expect_equal(colnames(ld3)[1], "F")
  #matrix Y with covar
  ld4 <- limma_dep(object=M, Y=y3, covar=covar)
  expect_equal(mean(ld3$p == ld4$p), 0)
})