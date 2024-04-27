des1 <- cbind(Int=1, design[, 1, drop=FALSE])
des2 <- cbind(Int=1, design[, 1, drop=FALSE], covar=1:6)

test_that("one phenotype", {
  # covar
  pheno.tab <- design[, 1, drop=FALSE]
  rc.des1 <- roast_cor(M, G=G, feat.tab=eztt, fun="fry", design = des1, prefix = "First3Arrays")
  rc.des2 <- roast_cor(M, G=G, feat.tab=eztt, fun="fry", design = des2, prefix = "First3Arrays")
  rmc.des1 <- roast_multi_cor(object=M, G=G, feat.tab=eztt, fun="fry", pheno.tab = pheno.tab)
  rmc.des2 <- roast_multi_cor(M, G=G, feat.tab=eztt, fun="fry", pheno.tab = pheno.tab, covariates = 1:6)
  expect_true(any(grepl("First3Arrays", colnames(rmc.des1))))
  expect_equal(colnames(rmc.des1), colnames(rmc.des2))
  expect_true(all(rmc.des1$p != rmc.des2$p))
  expect_equal(rc.des1, rmc.des1)
  expect_equal(rc.des2, rmc.des2)
  # weights
  res.rc <- roast_cor(M, G=G, feat.tab=eztt, pheno=pheno.v, fun="fry", weights=1:6, prefix = "abc")
  res.rmc <- roast_multi_cor(M, G=G, feat.tab=eztt, pheno.tab=data.frame(abc = pheno.v), fun="fry", weights=1:6, check.names = FALSE)
  expect_true(all(grepl(pattern = "abc", x = colnames(res.rmc)[-1])))
  expect_equal(res.rc, res.rmc)
  # NA
  pheno2.df <- data.frame(pheno2 = pheno2)
  expect_message(rcr.mw2 <- roast_cor(M, G=G, feat.tab=eztt, pheno=pheno2, fun="mroast", weights=1:6, prefix="pheno2"))
  expect_message(rmcr.mw2 <- roast_multi_cor(M, G=G, feat.tab=eztt, pheno.tab=pheno2.df, fun="mroast", weights=1:6))
  expect_equal(rcr.mw2, rmcr.mw2)
  # fun
  rmc.des3 <- roast_multi_cor(object=M, G=G, feat.tab=eztt, pheno.tab = pheno.tab)
  expect_equal(rmc.des3, rmc.des1)
})

test_that("multiple pheno", {
  # NA
  mph <- data.frame(design, pheno2=pheno2)
  # expect_message(rmc.multi <- roast_multi_cor(object=M, G=G, feat.tab=eztt, pheno.tab=mph, name = "mph")) # to test writing to file
  expect_message(rmc.multi <- roast_multi_cor(object=M, G=G, feat.tab=eztt, pheno.tab=mph, name = "mph"))
  rc.des1 <- roast_cor(M, G=G, feat.tab=eztt, fun="fry", design = des1, prefix = "First3Arrays")
  expect_equal(rmc.multi[, 1:6], rc.des1)
  expect_message(rc.pheno2 <- roast_cor(object=M, G=G, feat.tab=eztt, pheno=pheno2, prefix = "pheno2"))
  expect_equal(rmc.multi[, grep("pheno2", colnames(rmc.multi))], rc.pheno2[rownames(rmc.multi), -1])
  
  # covar & weights
  expect_message(rmc.multi2 <- roast_multi_cor(object=M, G=G, feat.tab=eztt, pheno.tab=mph, covariates = covar, weights=1:6))
  des.last3 <- model.matrix(~1 + design[,2] + covar)
  rc.des2b <- roast_cor(object=M, G=G, feat.tab=eztt, design = des.last3, prefix = "Last3Arrays", weights=1:6)
  expect_equal(rmc.multi2[, grep("Last3Arrays", colnames(rmc.multi2))], rc.des2b[rownames(rmc.multi2), -1])
  
  ## NA & covar & weights
  des.ph2 <- model.matrix(~1 + pheno2 + covar)
  ## rm first weight b/c pheno2[1] is NA
  rc.pheno2b <- roast_cor(object=M[, !is.na(pheno2)], G=G, feat.tab=eztt, design = des.ph2, prefix = "pheno2", weights=2:6)
  expect_equal(rmc.multi2[, grep("pheno2", colnames(rmc.multi2))], rc.pheno2b[rownames(rmc.multi2), -1])
})
