context("ezpcor")



test_that("ezpcor", {
  res1 <- ezcor(M, pheno.v, method="pearson", reorder.rows = FALSE)
  resp1 <- ezpreg(M, pheno.v, covar=numeric(length(pheno.v)), reorder.rows = FALSE)
  expect_equal(res1[,2:3], resp1[,2:3])
  
  ezpc <- limma_pcor(object=M, phenotype=pheno.v, covar=covar)
})