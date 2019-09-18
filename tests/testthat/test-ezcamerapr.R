context("ez cameraPR")

test_that("one comparison", {
  tmp1 <- ezcamerapr(eztt[,"Last3vsFirst3.logFC",drop=F], G=G, feat.tab=eztt, name = NA, adjust.method = "BH", 
                     alternative = "Up", min.nfeats = 3, max.nfeats = 100)
  expect_equal(tmp1[1,1], 10)
  expect_true(all.equal(tmp1[,2], sort(tmp1[,2], decreasing = TRUE)))
  expect_true(all.equal(order(tmp1[,3]), 1:nrow(tmp1)))
  expect_equal(rownames(tmp1)[3], "pwy1")
})

test_that("pwy1", {
  tmp2 <- ezcamerapr(eztt[,"Last3vsFirst3.logFC",drop=F], G=G, feat.tab=eztt, name = NA, adjust.method = "BH", 
                     alternative = "two.sided", min.nfeats = 3, max.nfeats = 20)
  expect_equal(rownames(tmp2)[1], "pwy1")
  # pwy1
  stat.v <- stats::setNames(eztt[, "Last3vsFirst3.logFC"], nm=rownames(eztt))
  vs <- limma::cameraPR(statistic = stat.v, index=G[[1]]$genes)
  expect_equal(vs$PValue, tmp2["pwy1", 3])

  tmp3 <- ezcamerapr(eztt[,c("Last3vsFirst3.logFC", "Last3vsFirst3.p")], G=G, feat.tab=eztt, name = NA, adjust.method = "BH", 
                     alternative = "Down", min.nfeats = 3, max.nfeats = 100)
  expect_equal(tmp3[1,1], 10)
  expect_equal(ncol(tmp3), 7)
  expect_equal(rownames(tmp3)[1], "pwy1")
  # pwy1 first in p-values down
  expect_equal(order(tmp3[, 6])[1], 1)
  
  tmp4 <- ezcamerapr(eztt[,c("Last3vsFirst3.logFC", "Last3vsFirst3.p")], G=G, feat.tab=eztt, name = NA, adjust.method = "BH", 
                     alternative = "two.sided", min.nfeats = 3, max.nfeats = 100)
  expect_equal(tmp4[1,1], 10)
  expect_equal(ncol(tmp4), 7)
  expect_equal(rownames(tmp4)[1], "pwy1")
  expect_equal(vs$PValue, tmp4["pwy1", 3])
})
