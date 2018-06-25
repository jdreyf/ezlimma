context("limma_dep")

# lc <- limma_cor(M, phenotype = pheno.v) in helper
test_that("numeric covariate", {
  ld <- limma_dep(M, Y = pheno.v)
  expect_equal(ld$p, lc$p)
  
  ld2 <- limma_dep(M, Y = pheno.v, covar=covar)
  expect_equal(mean(ld2$p == lc$p), 0)
  
  #pheno2 has an NA
  expect_error(limma_dep(M, Y = pheno2))
  
  expect_error(limma_dep(M, Y = rep(1, ncol(M))))
})

test_that("non-numeric y", {
  expect_error(limma_dep(M, Y = rep("a", ncol(M))))
  expect_error(limma_dep(M, Y = c(letters[1:(ncol(M)-1)], NA)))
  
  ld <- limma_dep(object=M, Y=grp)
  expect_equal(ld$p, eztt[rownames(ld), "Last3vsFirst3.p"])
  
  ld2 <- limma_dep(object=M, Y=grp, covar=covar)
  expect_equal(mean(ld2$p == eztt[rownames(ld2), "Last3vsFirst3.p"]), 0)
})

test_that("verbose", {
  expect_message(limma_dep(M, Y = pheno.v, covar=covar, verbose = TRUE))
  expect_message(limma_dep(object=M, Y=grp, verbose = TRUE))
})
