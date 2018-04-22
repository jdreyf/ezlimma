context("ezcor")

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
  
  #pheno2 has NAs
  res2 <- ezcor(M, pheno2, method="pearson", reorder.rows = TRUE, adjust.method = "BY")
  #cor.test but not cor handles NAs
  m.r2 <- apply(M, 1, function(v){
    cor.test(v, pheno2)$p.value
  })
  expect_equal(res2[rownames(M), "p"], m.r2)
  #reorder
  expect_false(rownames(res2)[1] == "gene1")
  expect_equal(colnames(res2)[3], "BY")
})