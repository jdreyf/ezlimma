context("limma cor")

test_that("matches topTable(eBayes(lmFit(M)))", {
  design <- model.matrix(~1+pheno.v)
  fit <- lmFit(M, design)
  fit2 <- eBayes(fit)
  toptab <- topTable(fit2, coef=2, num=Inf)
  # expect_equal(rownames(res2), rownames(toptab))
  expect_equal(lc[rownames(toptab), "p"], toptab$P.Value)
})

test_that("trend has effect", {
  lc2 <- limma_cor(M, phenotype = pheno.v, trend=TRUE)
  expect_equal(mean(lc2$p == lc[rownames(lc2), "p"]), 0)
})

test_that("weights", {
  lc.w <- limma_cor(M, phenotype = pheno.v, weights = 1:ncol(M))
  expect_equal(mean(lc.w$p==lc$p), 0)
  
  lc.gw <- limma_cor(M, phenotype = pheno.v, weights = 1:nrow(M))
  expect_equal(mean(lc.gw$p==lc$p), 0)
  
  #create EList object
  lc.el <- limma_cor(el, phenotype = pheno.v)
  expect_equal(mean(lc.el$p==lc$p), 0)
  
  #object$weights are being ignored
  expect_warning(lc.elw <- limma_cor(el, phenotype = pheno.v, weights = NULL))
  expect_equal(mean(lc.elw$p==lc$p), 1)
  
  # errors
  expect_error(lc.err <- limma_cor(M, phenotype=pheno.v, weights = 4))
  expect_error(lc.err <- limma_cor(M, phenotype=pheno.v, weights = "d"))
  expect_error(lc.err <- limma_cor(M, phenotype=pheno.v, weights = letters[1:2]))
  expect_error(lc.err <- limma_cor(M, phenotype=pheno.v, weights = letters[1:ncol(M)]))
  expect_error(lc.err <- limma_cor(M, phenotype=pheno.v, weights = c(1:(ncol(M)-1), NA)))
})

test_that("throw error for NAs in pheno", {
  expect_error(limma_cor(M, pheno2, reorder.rows = FALSE))
})

test_that("Handles NAs in M", {
  mm <- M
  mm[1, 1] <- NA
  expect_silent(lc.na <- limma_cor(mm, pheno.v, reorder.rows = FALSE))
  expect_true(lc.na[1, "p"] != lc["gene1", "p"])
  
  mm[1, 1:5] <- NA
  expect_warning(lc.na <- limma_cor(mm, pheno.v, reorder.rows = TRUE))
  expect_true(is.na(lc.na["gene1", "p"]))
  
  mm[1, 1:6] <- NA
  expect_silent(lc.na <- limma_cor(mm, pheno.v, reorder.rows = TRUE))
  expect_true(is.na(lc.na["gene1", "p"]))
})

test_that("matches limma_contrast grp-means parametrization", {
  lc <- limma_cor(M, phenotype = design[,"Last3Arrays"])
  expect_equal(lc$p, eztt[rownames(lc), "Last3vsFirst3.p"])
})

test_that("matches limma_contrast with covariate", {
  des2 <- model.matrix(~0+grp+covar)
  colnames(des2) <- gsub("grp", "", colnames(des2))
  eztt2 <- limma_contrasts(object=M, design = des2, contrast.v = contr.v[3])
  
  des.lc <- model.matrix(~1+design[,"Last3Arrays"]+covar)
  lc2 <- limma_cor(M, design = des.lc)
  expect_equal(lc2$p, eztt2[rownames(lc2), "Last3vsFirst3.p"])
  
  #dupcor has effect
  eztt3 <- limma_contrasts(object=M, design = des2, contrast.v = contr.v[3], block = grp, correlation = 0.9)
  expect_lte(mean(eztt2$Last3vsFirst3.p == eztt3[rownames(eztt2), "Last3vsFirst3.p"]), 0.01)
  lc3 <- limma_cor(M, design = des.lc, block = grp, correlation = 0.9)
  expect_lte(mean(lc2$p == lc3[rownames(lc2), "p"]), 0.01)
  
  expect_equal(lc3$p, eztt3[rownames(lc3), "Last3vsFirst3.p"])
})

test_that("z col", {
  lc2 <- limma_cor(M, phenotype = pheno.v, trend=TRUE, cols = "z")
  expect_equal(colnames(lc2), "z")
  
  lc3 <- limma_cor(M, phenotype = pheno.v, trend=TRUE, cols = c("z", "P.Value"))
  expect_equal(colnames(lc3)[1], "z")
  expect_equal(colnames(lc3)[2], "p")
})

test_that("!moderated", {
  expect_error(limma_cor(object=M, phenotype = pheno.v, trend=TRUE, moderated = FALSE))
  
  lc2 <- limma_cor(object=M, phenotype = pheno.v, moderated = FALSE,
                   cols=c("logFC", "t", "P.Value", "adj.P.Val"))
  g1.t <- cor.test(x=M["gene1", ], y=pheno.v)
  expect_equal(lc2["gene1", "t"], unname(g1.t$statistic))
  expect_equal(lc2["gene1", "p"], g1.t$p.value)
  
  # with covariates
  des.lc <- model.matrix(~1+pheno.v+covar)
  lc <- limma_cor(M[1,, drop=FALSE], design = des.lc)
  lc2 <- limma_cor(object=M[1,, drop=FALSE], design = des.lc, moderated = FALSE)
  lc3 <- limma_cor(object=M, design = des.lc, moderated = FALSE)
  expect_equal(lc$p, lc2$p)
  expect_equal(lc$p, lc3["gene1", "p"])
  
  # independent of which of M or Y is on LHS
  # so can use limma_cor without concern for which var to put on LHS!
  des.lc2 <- model.matrix(~1+M[1,]+covar)
  pheno.mat <- matrix(pheno.v, ncol=length(pheno.v), dimnames=list("pheno", names(pheno.v)))
  lc4 <- limma_cor(object=pheno.mat, design = des.lc2, moderated = FALSE)
  expect_equal(lc2$p, lc4$p)
  
  # matches ppcor
  # pc <- ppcor::pcor.test(x=pheno.v, y=M[1,], z=covar)
  # pc2 <- ppcor::pcor.test(y=pheno.v, x=M[1,], z=covar)
  # expect_equal(pc$p.value, pc2$p.value)
  # pc$p.value is 0.9643012
  expect_equal(signif(lc4[1, "p"], digits = 4), 0.9643)
})

test_that("reduce.df has effect", {
  lc2 <- limma_cor(M, phenotype = pheno.v, reduce.df = 3)
  expect_gte(mean(lc2$p >= lc[rownames(lc2), "p"]), 0.6)
})

test_that("size", {
  lc1 <- limma_cor(M[-1,], phenotype = covar)
  lc2 <- limma_cor(M[-1,], phenotype = phenotype)
  set.seed(42)
  ph3 <- rnorm(ncol(M))
  lc3 <- limma_cor(M[-1,], phenotype = ph3)
  ph4 <- rnorm(ncol(M))
  lc4 <- limma_cor(M[-1,], phenotype = ph4)
  
  expect_lte(mean(rbind(lc1, lc2, lc3, lc4)[,"p"] < 0.05), 0.05)
})

test_that("ndups & spacing", {
  expect_error(limma_cor(M, phenotype = pheno.v, block=rep(1:3, times=2), ndups=2))
  # need cor if ndups=2 or !is.null(block) 
  expect_error(limma_cor(M, phenotype = pheno.v, ndups=2))
  
  lc.nodups <- limma_cor(object=M, phenotype = pheno.v, ndups = 1, spacing=1)
  expect_true(all.equal(lc, lc.nodups))
  
  ph2 <- phenotype
  ph2[1:3] <- ph2[1:3]-4
  
  lc.dups <- limma_cor(object=M, phenotype = ph2, ndups = 2, correlation = 0.3)
  lc.dups2 <- limma_cor(object=M, phenotype = ph2, ndups = 2, correlation = 0.6)
  expect_equal(nrow(lc.dups), 50)
  expect_equal(dim(lc.dups), dim(lc.dups2))
  expect_true(!all(lc.dups == lc.dups2))
  expect_equal(nrow(lc.dups), 50)
  
  # gene numbers in lc.dups are all odds from 1:99:2
  expect_equal(sort(as.numeric(sub("gene", "", rownames(lc.dups)))), seq(1, 99, by=2))
  expect_equal(lc.dups["gene1", "AveExpr"], mean(M[1:2, 1:6]))
  expect_gt(lc.dups["gene1", "AveExpr"], 0)
  expect_lt(abs(mean(lc.dups[, "AveExpr"])), 0.1)
  
  expect_lt(lc.dups["gene1", "slope"], 0)
  expect_lt(lc.dups2["gene1", "slope"], 0)
})