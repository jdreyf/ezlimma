context("modify_hitman_pvalues")

hm <- hitman(E=grp, M=M, Y=M[1,])
hm$MY2.p <- hm$MY.p
hm$EM2.p <- hm$EM.p
ey.sign <- sign(limma_dep(object=M[1,], Y=grp, covariates=NULL, prefix="EY")$EY.t)
p.cols <- c("EM2.p", "MY2.p")
stat.cols <- c("EM.t", "MY.slope")
ret <- modify_hitman_pvalues(tab=hm, overall.sign = ey.sign, p.cols=p.cols)

#gene at bottom of list
#my.p > em.p
test_that("gene88", {
  v <- ret["gene88",] #a matrix
  prod.sign <- sign(v[,stat.cols[1]])*sign(v[,stat.cols[2]])
  expect_false(prod.sign == ey.sign)
  
  expect_gt(v[,"MY.slope"], 0)
  expect_gt(v[,"MY.p"], v[,"EM.p"])
  expect_equal(v[,"EM.p"], v[,"EM2.p"])
  alt.tmp <- ifelse (v[,stat.cols[1]] > 0, "less", "greater")
  expect_equal(as.numeric(v["MY2.p"]), as.numeric(two2one_tailed(v, stat.col = "MY.slope", p.col="MY.p", alternative = alt.tmp)))
})

#my.p < em.p, & wrong sign
test_that("gene72", {
  v <- ret["gene72",] #a matrix
  prod.sign <- sign(v[,stat.cols[1]])*sign(v[,stat.cols[2]])
  expect_false(prod.sign == ey.sign)
  
  expect_lt(v[,"MY.slope"], 0)
  expect_lt(v[,"MY.p"], v[,"EM.p"])
  expect_equal(v[,"MY.p"], v[,"MY2.p"])
  alt.tmp <- ifelse (v[,stat.cols[1]] > 0, "less", "greater")
  expect_equal(as.numeric(v["EM2.p"]), as.numeric(two2one_tailed(v, stat.col = "EM.t", p.col="EM.p", alternative = alt.tmp)))
})