context("ezpcor")

test_that("no covar", {
  res0 <- limma_pcor(object=M, phenotype=pheno.v, batch=grp)
  res1 <- limma_pcor(object=M, phenotype=pheno.v, covariates = as.numeric(as.factor(grp)))
  # expect_equal(res0, res1)
  
  m2 <- rbind(M, phenotype)
  res2 <- limma_pcor(object=m2, phenotype=pheno.v, batch=grp)
  expect_equal(rownames(res2)[1], "phenotype")
})

test_that("with covar", {
  res0 <- limma_pcor(object=M, phenotype=pheno.v, batch=grp, covariates = covar)
  res1 <- limma_pcor(object=M, phenotype=pheno.v, batch=as.numeric(as.factor(grp)), covariates = covar)
  expect_equal(res0, res1)
  
  m2 <- rbind(M, phenotype, covariates = covar)
  res2 <- limma_pcor(object=m2, phenotype=pheno.v, batch=grp, covariates = covar)
  #possible that phenotype wouldn't be first, but it is most likely
  expect_equal(rownames(res2)[1], "phenotype")
  
  expect_error(res2 <- limma_pcor(object=m2, phenotype=pheno.v, batch=grp, covariates = grp))
})
