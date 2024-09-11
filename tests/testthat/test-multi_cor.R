context("multi_cor")

test_that("matches ezcor", {
  res.ez <- ezcor(M, pheno2, method="spearman", reorder.rows = TRUE)
  expect_lte(res.ez[1, "p"], min(res.ez[-1, "p"]))
  res.mc <- multi_cor(M, cbind(a=pheno.v, b=pheno2), method="spearman", reorder.rows = FALSE)
  expect_equal(as.numeric(res.ez[rownames(M), "p"]), as.numeric(res.mc[,"b.p"]))
  expect_error(multi_cor(M, pheno.v, method="limma", reorder.rows = FALSE))
})

test_that("matches limma", {
  res.lc <- limma_cor(M, pheno.v, reorder.rows = FALSE, prefix = colnames(pheno.v), cols = c('AveExpr', 'P.Value', 'adj.P.Val', 'logFC'))
  res.mc2 <- multi_cor(M, cbind(a=pheno.v, b=pheno.v), method="limma", reorder.rows = FALSE)
  res.lc <- res.lc[rownames(res.mc2),]
  res.lc.pval <- setNames(res.lc$p,rownames(res.lc))
  expect_equal(as.numeric(res.lc.pval), as.numeric(res.mc2[,"b.p"]))
  
  # with NAs
  pheno2.nona <- na.omit(pheno2)
  expect_error(res.mc.na <- multi_cor(M, cbind(a=pheno2, b=pheno.v, c=rep(NA, length(pheno.v))), method="limma", reorder.rows = FALSE))
  res.lc2 <- limma_cor(M[, !is.na(pheno2)], pheno2.nona, reorder.rows = FALSE, prefix = colnames(pheno.v), cols = c('AveExpr', 'P.Value', 'adj.P.Val', 'logFC'))
  res.mc3 <- multi_cor(M, cbind(a=pheno2, b=pheno.v), method="limma", reorder.rows = FALSE)
  expect_equal(res.mc3[rownames(res.mc2), "b.p"], res.mc2$b.p)
  expect_equal(res.mc3[rownames(res.lc2), "a.p"], res.lc2$p)
  
  # with block
  block.tmp <- c("a", letters[1:5])
  res.lc3 <- limma_cor(M, phenotype = pheno.v, block = block.tmp, correlation = -0.5, reorder.rows = FALSE)
  res.mc4 <- multi_cor(M, cbind(a=pheno2, b=pheno2, c=pheno.v, d=pheno.v), method="limma", block = block.tmp, 
                       correlation = c(0, 0.9, -0.5, 0.5), reorder.rows = FALSE)
  expect_equal(res.lc3$p, res.mc4$c.p)
  expect_equal(res.lc3$slope, res.mc4$c.slope)
  # block should matter for pheno.v but not for pheno2, since it only one sample per block
  expect_equal(res.mc4$a.p, res.mc4$b.p)
  expect_true(res.mc4$c.p[1] != res.mc4$d.p[1])
})

test_that("rows get reordered", {
  res.mc2 <- multi_cor(M, cbind(a=pheno.v, b=pheno.v), method="limma", reorder.rows = FALSE)
  res.mc3 <- multi_cor(M, cbind(a=pheno.v, b=pheno.v), method="limma", reorder.rows = TRUE)
  expect_false(all(order(rownames(res.mc3)) == order(rownames(res.mc2))))
})

test_that("covars", {
  des.lc <- model.matrix(~1+design[,"Last3Arrays"]+covar)
  lc <- limma_cor(M, design = des.lc, prefix = "Last3Arrays")
  mc <- multi_cor(object=M, pheno.tab = design[,"Last3Arrays", drop=FALSE], method = "limma", covariates = covar)
  expect_equal(mc, lc)
  
  covar.mat <- cbind(covar, rnorm(length(covar)))
  des.lc <- model.matrix(~1+design[,"Last3Arrays"]+covar.mat)
  lc <- limma_cor(M, design = des.lc, prefix = "Last3Arrays")
  mc <- multi_cor(object=M, pheno.tab = design[,"Last3Arrays", drop=FALSE], method = "limma", covariates = covar.mat)
  expect_equal(mc, lc)
  
  covar.df <- data.frame(covar, rnorm(length(covar)))
  expect_error(multi_cor(object=M, pheno.tab = design[,"Last3Arrays", drop=FALSE], method = "limma", covariates = covar.df))
  
  expect_warning(multi_cor(object=M, pheno.tab = design[,"Last3Arrays", drop=FALSE], covariates = covar))
  
  # w/ NAs
  des.lc2 <- model.matrix(~1+pheno2+covar)
  lc2 <- limma_cor(M[, -1], design = des.lc2, prefix = "pheno2")
  mc2 <- multi_cor(object=M, pheno.tab = data.frame(pheno2), method = "limma", covariates = covar, check.names = FALSE)
  expect_equal(mc2, lc2)
})