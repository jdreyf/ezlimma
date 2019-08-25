context("ez cameraPR")

test_that("alternative", {
  tmp1 <- ezcamerapr(gstats=eztt[,"Last3vsFirst3.logFC",drop=F], G=G, feat.tab=eztt, name = NA, adjust.method = "BH", 
                     alternative = "Up", min.nfeats = 3, max.nfeats = 100)
  expect_equal(tmp1$NGenes[1], 10)
  expect_true(all.equal(tmp1$Direction, sort(tmp1$Direction, decreasing = TRUE)))
  expect_true(all.equal(order(tmp1$p), 1:nrow(tmp1)))
  expect_equal(rownames(tmp1)[3], "pwy1")
})

test_that("pwy1", {
  tmp2 <- ezcamerapr(gstats=eztt[,"Last3vsFirst3.logFC",drop=F], G=G, feat.tab=eztt, name = NA, adjust.method = "BH", 
                     alternative = "two.sided", min.nfeats = 3, max.nfeats = 20)
  expect_equal(rownames(tmp2)[1], "pwy1")
  # pwy1
  stat.v <- stats::setNames(eztt[, "Last3vsFirst3.logFC"], nm=rownames(eztt))
  vs <- limma::cameraPR(statistic = stat.v, index=G[[1]]$genes)
  expect_equal(vs$PValue, tmp2["pwy1", "p"])
})