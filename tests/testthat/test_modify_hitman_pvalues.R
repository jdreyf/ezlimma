context("modify_hitman_pvalues")

hm <- hitman(E=grp2, M=M, Y=M[1,])
hm$MY_dir.p <- hm$MY.p
hm$EM_dir.p <- hm$EM.p
ey.sign <- sign(limma_dep(object=M[1,], Y=grp2, covariates=NULL, prefix="EY")$EY.t)
p.cols <- c("EM_dir.p", "MY_dir.p")
ret <- modify_hitman_pvalues(tab=hm, overall.sign = ey.sign, p.cols=p.cols)

#gene at bottom of list
#my.p > em.p & wrong sign --> flip my.p sign
test_that("gene88", {
  v <- ret["gene88",] #a matrix
  prod.sign <- sign(v[,"EM.t"])*sign(v[,"MY.slope"])
  expect_false(prod.sign == ey.sign)
  
  expect_gt(v[,"MY.slope"], 0)
  expect_gt(v[,"MY.p"], v[,"EM.p"])
  expect_equal(v[,"EM.p"], v[,"EM_dir.p"])
  alt.tmp <- ifelse (v[,"EM.t"] > 0, "greater", "less")
  expect_equal(as.numeric(v["MY_dir.p"]), 
               as.numeric(two2one_tailed(v, stat.col = "MY.slope", p.col="MY.p", alternative = alt.tmp)))
})

#my.p < em.p, & wrong sign
test_that("gene72", {
  v <- ret["gene72",] #a matrix
  prod.sign <- sign(v[,"EM.t"])*sign(v[,"MY.slope"])
  expect_false(prod.sign == ey.sign)
  
  expect_lt(v[,"MY.slope"], 0)
  expect_lt(v[,"MY.p"], v[,"EM.p"])
  expect_equal(v[,"MY.p"], v[,"MY_dir.p"])

  alt.tmp <- ifelse (v[,"EM.t"] > 0, "less", "greater")
  expect_equal(as.numeric(v["EM_dir.p"]), 
               as.numeric(two2one_tailed(v, stat.col = "EM.t", p.col="EM.p", alternative = alt.tmp)))
})