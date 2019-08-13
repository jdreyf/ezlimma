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

test_that("z col", {
  lc2 <- limma_cor(M, phenotype = pheno.v, trend=TRUE, cols = "z")
  expect_equal(colnames(lc2), "z")
  
  lc3 <- limma_cor(M, phenotype = pheno.v, trend=TRUE, cols = c("z", "P.Value"))
  expect_equal(colnames(lc3)[1], "z")
  expect_equal(colnames(lc3)[2], "p")
})

test_that("!moderated", {
  expect_error(limma_cor(object=M, phenotype = pheno.v, trend=TRUE, moderated = FALSE))
  
  lc2 <- limma_cor(object=M, phenotype = pheno.v, moderated = FALSE,
                   cols=c("logFC", "t", "P.Value", "adj.P.Val"))
  g1.t <- cor.test(x=M["gene1", ], y=pheno.v)
  expect_equal(lc2["gene1", "t"], unname(g1.t$statistic))
  expect_equal(lc2["gene1", "p"], g1.t$p.value)
  
  # with covariates
  des.lc <- model.matrix(~1+pheno.v+covar)
  lc <- limma_cor(M[1,], design = des.lc)
  lc2 <- limma_cor(object=M[1,], design = des.lc, moderated = FALSE)
  lc3 <- limma_cor(object=M, design = des.lc, moderated = FALSE)
  expect_equal(lc$p, lc2$p)
  expect_equal(lc$p, lc3["gene1", "p"])
  
  # independent of which of M or Y is on LHS
  # so can use limma_cor without concern for which var to put on LHS!
  des.lc2 <- model.matrix(~1+M[1,]+covar)
  lc4 <- limma_cor(object=pheno.v, design = des.lc2, moderated = FALSE)
  expect_equal(lc2$p, lc4$p)
  
  # matches ppcor
  pc <- ppcor::pcor.test(x=pheno.v, y=M[1,], z=covar)
  pc2 <- ppcor::pcor.test(y=pheno.v, x=M[1,], z=covar)
  expect_equal(pc$p.value, pc2$p.value)
  expect_equal(lc4[1, "p"], pc$p.value)
})

test_that("reduce.df has effect", {
  lc2 <- limma_cor(M, phenotype = pheno.v, reduce.df = 3)
  expect_gte(mean(lc2$p >= lc[rownames(lc2), "p"]), 0.6)
})

test_that("size", {
  lc1 <- limma_cor(M[-1,], phenotype = covar)
  lc2 <- limma_cor(M[-1,], phenotype = phenotype)
  set.seed(42)
  ph3 <- rnorm(ncol(M))
  lc3 <- limma_cor(M[-1,], phenotype = ph3)
  ph4 <- rnorm(ncol(M))
  lc4 <- limma_cor(M[-1,], phenotype = ph4)
  
  expect_lte(mean(rbind(lc1, lc2, lc3, lc4)[,"p"] < 0.05), 0.05)
})