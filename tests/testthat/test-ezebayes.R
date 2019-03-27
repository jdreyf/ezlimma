context("ezebayes")

test_that("ezebayes", {
  object=M; contrast.v = contr.v
  design=NULL;weights=NA;block=NULL;
  correlation=NULL;adjust.method="BH";add.means=!is.null(grp);treat.lfc=NULL;
  cols=c("P.Value","adj.P.Val","logFC")
  design <- stats::model.matrix(~0+grp)
  colnames(design) <- sub("grp", "", colnames(design), fixed=TRUE)
  fit <- limma::lmFit(object, design=design, block = block, correlation = correlation)
  contr.mat <- limma::makeContrasts(contrasts=contrast.v, levels=design)
  fit <- limma::contrasts.fit(fit, contr.mat)
  fit2.ez <- ezebayes(fit = fit, moderated = FALSE)
  
  # warning when !moderated & trend
  expect_warning(ezebayes(fit = fit, moderated = FALSE, trend = TRUE))
  
  # get correct components
  expect_true("t" %in% names(fit2.ez))
  expect_true("F" %in% names(fit2.ez))
  
  # get values as calculated
  g1.lmfit.t <- (fit$coefficients / fit$stdev.unscaled / fit$sigma)[1, 3]
  expect_equal(fit2.ez$t["gene1", 3], g1.lmfit.t)
  
  # different from eBayes when !moderated
  fit2.eb <- ezebayes(fit = fit, moderated = TRUE)
  expect_true(all(fit2.eb$t != fit2.ez$t))
  expect_true(all(fit2.eb$F != fit2.ez$F))
  expect_false("df.total" %in% names(fit2.ez))
  
  # t & p = OLS estimate
  # only test two-group comparison, since limma looks at all samples to cal var
  g1.t <- t.test(x=M["gene1", grp=="Last3"], y=M["gene1", grp=="First3"], var.equal = TRUE)
  expect_equal(setNames(fit2.ez$t["gene1", 3], nm="t"), g1.t$statistic)
})