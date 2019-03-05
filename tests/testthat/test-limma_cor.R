context("limma cor")

test_that("matches topTable(eBayes(lmFit(M)))", {
  design <- model.matrix(~1+pheno.v)
  fit <- lmFit(M, design)
  fit2 <- eBayes(fit)
  toptab <- topTable(fit2, coef=2, num=Inf)
  # expect_equal(rownames(res2), rownames(toptab))
  expect_equal(lc[rownames(toptab), "p"], toptab$P.Value)
})

test_that("trend has effect", {
  lc2 <- limma_cor(M, phenotype = pheno.v, trend=TRUE)
  expect_equal(mean(lc2$p == lc[rownames(lc2), "p"]), 0)
})

test_that("weights", {
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

test_that("throw error for NAs in pheno", {
  expect_error(limma_cor(M, pheno2, reorder.rows = FALSE))
})

test_that("matches limma_contrast grp-means parametrization", {
  lc <- limma_cor(M, phenotype = design[,"Last3Arrays"])
  expect_equal(lc$p, eztt[rownames(lc), "Last3vsFirst3.p"])
})

test_that("matches limma_contrast with covariate", {
  des2 <- model.matrix(~0+grp+covar)
  colnames(des2) <- gsub("grp", "", colnames(des2))
  eztt2 <- limma_contrasts(object=M, design = des2, contrast.v = contr.v[3])
  des.lc <- model.matrix(~1+design[,"Last3Arrays"]+covar)
  lc2 <- limma_cor(M, design = des.lc)
  expect_equal(lc2$p, eztt2[rownames(lc2), "Last3vsFirst3.p"])
  
  #dupcor has effect
  eztt3 <- limma_contrasts(object=M, design = des2, contrast.v = contr.v[3], block = grp, correlation = 0.9)
  expect_lte(mean(eztt2$Last3vsFirst3.p == eztt3[rownames(eztt2), "Last3vsFirst3.p"]), 0.01)
  lc3 <- limma_cor(M, design = des.lc, block = grp, correlation = 0.9)
  expect_lte(mean(lc2$p == lc3[rownames(lc2), "p"]), 0.01)
  
  expect_equal(lc3$p, eztt3[rownames(lc3), "Last3vsFirst3.p"])
})