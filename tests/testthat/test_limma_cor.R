context("limma cor")

test_that("limma_cor matches topTable(eBayes(lmFit(M)))", {
  design <- model.matrix(~1+pheno.v)
  fit <- lmFit(M, design)
  fit2 <- eBayes(fit)
  toptab <- topTable(fit2, coef=2, num=Inf)
  # expect_equal(rownames(res2), rownames(toptab))
  expect_equal(lc[rownames(toptab), "p"], toptab$P.Value)
})

test_that("limma_cor trend has effect", {
  lc2 <- limma_cor(M, phenotype = pheno.v, trend=TRUE)
  expect_equal(mean(lc2$p == lc[rownames(lc2), "p"]), 0)
})

test_that("limma_cor weights", {
  lc.w <- limma_cor(M, phenotype = pheno.v, weights = 1:ncol(M))
  expect_equal(mean(lc.w$p==lc$p), 0)
  
  lc.gw <- limma_cor(M, phenotype = pheno.v, weights = 1:nrow(M))
  expect_equal(mean(lc.gw$p==lc$p), 0)
  
  #create EList object
  lc.el <- limma_cor(el, phenotype = pheno.v)
  expect_equal(mean(lc.el$p==lc$p), 0)
  
  #object$weights are being ignored
  expect_warning(lc.elw <- limma_cor(el, phenotype = pheno.v, weights = NULL))
  expect_equal(mean(lc.elw$p==lc$p), 1)
})
