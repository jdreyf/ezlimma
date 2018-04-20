context("cor")

test_that("ezcor with different methods matches cor", {
  #make separate apply calls to ensure no bugs
  m.r <- apply(M, 1, function(v){
    cor(v, pheno.v)
  })
  m.p <- apply(M, 1, function(v){
    cor.test(v, pheno.v)$p.value
  })
  res1 <- ezcor(M, pheno.v, method="pearson", reorder.rows = FALSE)
  expect_equal(res1[,"cor"], m.r)
  expect_equal(res1[,"p"], m.p)
  
  #make separate apply calls to ensure no bugs
  m.r <- apply(M, 1, function(v){
    cor(v, pheno.v, method = "kendall")
  })
  m.p <- apply(M, 1, function(v){
    cor.test(v, pheno.v, method = "kendall")$p.value
  })
  res1 <- ezcor(M, pheno.v, method="kendall", reorder.rows = FALSE)
  expect_equal(res1[,"tau"], m.r)
  expect_equal(res1[,"p"], m.p)
})

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

test_that("multi_cor matches ezcor & limma_cor", {
  res.ez <- ezcor(M, pheno2, method="spearman", reorder.rows = FALSE)
  res.mc <- multi_cor(M, cbind(a=pheno.v, b=pheno2), method="spearman", reorder.rows = FALSE)
  expect_equal(res.ez[,"p"], res.mc[,"b.p"])
  
  res.lm <- data.matrix(limma_cor(M, pheno2, reorder.rows = FALSE))
  res.mc <- multi_cor(M, cbind(a=pheno.v, b=pheno2), method="limma", reorder.rows = FALSE)
  expect_equal(res.lm[,"p"], res.mc[,"b.p"])
})
