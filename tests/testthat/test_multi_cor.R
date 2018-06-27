context("multi_cor")

test_that("matches ezcor", {
  res.ez <- ezcor(M, pheno2, method="spearman", reorder.rows = TRUE)
  expect_lte(res.ez[1, "p"], min(res.ez[-1, "p"]))
  res.mc <- multi_cor(M, cbind(a=pheno.v, b=pheno2), method="spearman", reorder.rows = FALSE)
  expect_equal(res.ez[rownames(M), "p"], res.mc[,"b.p"])
  expect_error(multi_cor(M, pheno.v, method="limma", reorder.rows = FALSE))
})

test_that("matches limma", {
  res.lc <- limma_cor(M,pheno.v, reorder.rows = FALSE,prefix = colnames(pheno.v), cols = c('AveExpr', 'P.Value', 'adj.P.Val', 'logFC'))
  res.mc2 <- multi_cor(M, cbind(a=pheno.v, b=pheno.v), method="limma", reorder.rows = FALSE)
  res.lc <- res.lc[rownames(res.mc2),]
  res.lc.pval <- setNames(res.lc$p,rownames(res.lc))
  expect_equal(res.lc.pval, res.mc2[,"b.p"])
})

test_that("rows get reordered", {
  res.mc2 <- multi_cor(M, cbind(a=pheno.v, b=pheno.v), method="limma", reorder.rows = FALSE)
  res.mc3 <- multi_cor(M, cbind(a=pheno.v, b=pheno.v), method="limma", reorder.rows = TRUE)
  expect_false(all(order(rownames(res.mc3)) == order(rownames(res.mc2))))
})