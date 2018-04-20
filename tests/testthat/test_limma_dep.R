context("limma_dep")


# lc <- limma_cor(M, phenotype = pheno.v) in helper
test_that("limma_dep with numeric covariate", {
  ld <- limma_dep(M, y = pheno.v)
  expect_equal(ld$p, lc$p)
  
  ld2 <- limma_dep(M, y = pheno.v, covar=covar)
  expect_equal(mean(ld2$p == lc$p), 0)
  
  #pheno2 has an NA
  expect_message(lc2 <- limma_cor(M, phenotype = pheno2))
  ld3 <- limma_dep(M, y = pheno2)
  expect_equal(ld3$p, lc2$p)
  
  expect_error(limma_dep(M, y = rep(1, ncol(M))))
})

test_that("limma_dep with non-numeric y", {
  expect_error(limma_dep(M, y = rep("a", ncol(M))))
  expect_error(limma_dep(M, y = c(letters[1:(ncol(M)-1)], NA)))
  
  ld <- limma_dep(object=M, y=grp)
  expect_equal(ld$p, eztt[rownames(ld), "Last3vsFirst3.p"])
  
  ld2 <- limma_dep(object=M, y=grp, covar=covar)
  expect_equal(mean(ld2$p == eztt[rownames(ld2), "Last3vsFirst3.p"]), 0)
})
