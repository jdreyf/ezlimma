context("limmaF")

des2 <- cbind(intercept=rep(1,6), Last3Arrays=c(0,0,0,1,1,1), covar=covar)

test_that("limmaF", {
  expect_message(lf <- limmaF(object=M, design=des2[,1:2]))
  
  lf2 <- limmaF(object=M, design=des2[,1:2], coef=colnames(des2)[2])
  expect_equal(lf2$p, eztt[rownames(lf2), "Last3vsFirst3.p"])
  expect_equal(rownames(lf2)[1], "gene1")
  
  lf3 <- limmaF(object=M, design=des2, coef=colnames(des2)[2:3])
  expect_equal(rownames(lf3)[1], "gene1")
  expect_equal(mean(lf2$p==lf3$p), 0)
  
  lf4 <- limmaF(object=M, design=des2[,1:2], coef=colnames(des2)[2], weights=runif(ncol(M)), reorder.rows = FALSE, prefix = "wts")
  expect_equal(rownames(lf4)[2], "gene2")
  expect_equal(colnames(lf4), c("wts.t", "wts.p"))
  expect_equal(mean(lf2$p==lf4$wts.p), 0)
})