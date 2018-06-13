context("multi_cor")

test_that("matches ezcor", {
  res.ez <- ezcor(M, pheno2, method="spearman", reorder.rows = TRUE)
  expect_lte(res.ez[1, "p"], min(res.ez[-1, "p"]))
  res.mc <- multi_cor(M, cbind(a=pheno.v, b=pheno2), method="spearman", reorder.rows = FALSE)
  expect_equal(res.ez[rownames(M), "p"], res.mc[,"b.p"])
  
  expect_error(multi_cor(M, pheno.v, method="limma", reorder.rows = FALSE))
})
