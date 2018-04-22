context("multi_cor")

test_that("multi_cor matches ezcor & limma_cor", {
  res.ez <- ezcor(M, pheno2, method="spearman", reorder.rows = TRUE)
  expect_lte(res.ez[1, "p"], min(res.ez[-1, "p"]))
  res.mc <- multi_cor(M, cbind(a=pheno.v, b=pheno2), method="spearman", reorder.rows = FALSE)
  expect_equal(res.ez[rownames(M), "p"], res.mc[,"b.p"])
  
  res.lm <- data.matrix(limma_cor(M, pheno2, reorder.rows = FALSE))
  res.mc <- multi_cor(M, cbind(a=pheno.v, b=pheno2), method="limma", reorder.rows = FALSE)
  expect_equal(res.lm[,"p"], res.mc[,"b.p"])
})
