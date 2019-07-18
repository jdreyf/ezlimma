context("roast contrasts")

test_that("weights", {
  #aw have effect
  rcn.fw <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="fry", weights=1:6)
  expect_equal(mean(rcn.f$First3.p==rcn.fw$First3.p), 0)
  rcn.mw <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast", weights=1:6)
  expect_equal(mean(rcn.m$First3.p==rcn.mw$First3.p), 0)
  
  #gene weights have effect
  expect_warning(rcn.fw <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="fry", 
                                           gene.weights=1:nrow(M)))
  expect_equal(mean(rcn.f$First3.p==rcn.fw$First3.p), 1)
  rcn.mw <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast", gene.weights=1:nrow(M))
  expect_equal(mean(rcn.m$First3.p==rcn.mw$First3.p), 0)
  
  rcn.fe <- roast_contrasts(el, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="fry")
  expect_equal(mean(rcn.f$First3.p==rcn.fe$First3.p), 0)
  rcn.me <- roast_contrasts(el, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast")
  expect_equal(mean(rcn.m$First3.p==rcn.me$First3.p), 0)
  
  #suppress object$weights
  expect_warning(rcn.fw <- roast_contrasts(el, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="fry", 
                                           weights=el$weights))
  expect_equal(mean(rcn.fw$First3.p==rcn.fe$First3.p), 1)
  expect_warning(rcn.mw <- roast_contrasts(el, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast", 
                                           weights=el$weights))
  expect_equal(mean(rcn.mw$First3.p==rcn.me$First3.p), 1)
})

test_that("one sided testing", {
  tmp <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="fry")
  tmp2 <- roast_contrasts(object=M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="fry", alternative = "less")
  expect_equal(1-tmp["pwy1", "First3.p"]/2, tmp2["pwy1", "First3.p"])
  expect_equal(tmp["pwy2", "First3.p"]/2, tmp2["pwy2", "First3.p"])
  expect_equal(1-tmp["pwy3", "First3.p"]/2, tmp2["pwy3", "First3.p"])
  #no mixed columns
  expect_equal(length(grep("Mixed", colnames(tmp2))), 0)
  
  tmp3 <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast")
  tmp4 <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast", alt="less")
  expect_lt(abs(1-tmp3["pwy1", "First3.p"]/2 - tmp4["pwy1", "First3.p"]), 0.001)
  expect_lt(abs(tmp3["pwy2", "First3.p"]/2 - tmp4["pwy2", "First3.p"]), 0.001)
  expect_lt(abs(1-tmp3["pwy3", "First3.p"]/2 - tmp4["pwy3", "First3.p"]), 0.001)
})

test_that("trend has effect", {
  rc.res1 <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="fry", weights=1:6)
  rc.res1t <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="fry", weights=1:6, trend = TRUE)
  expect_equal(mean(rc.res1$First3.p == rc.res1t$First3.p), 0)
  
  rc.res2 <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast", weights=1:6)
  rc.res2t <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast", weights=1:6, 
                              trend = TRUE)
  p.cols <- grep("\\.p", colnames(rc.res2))
  expect_equal(mean(rc.res2[,p.cols] == rc.res2t[,p.cols]), 0)
})

test_that("other values of adjust.method", {
  #permutations are comparable between adjust.methods b/c seed
  bonferroni  <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast", weights=1:6, 
                                 adjust.method = 'bonferroni')
  default <- roast_contrasts(M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast", weights=1:6)
  expect_false(all(colnames(default) == colnames(bonferroni)))
  
  expect_true(all(bonferroni$First3.bonferroni >= default[rownames(bonferroni), "First3.FDR"]))
})

test_that("proportion up & down", {
  #First3.PropDownP05 is all zero
  g.genes <- c(G[[1]]$genes, G[[2]]$genes, G[[3]]$genes)
  expect_equal(sum(eztt[g.genes, "First3.p"] <= 0.05 & eztt[g.genes, "First3.logFC"] < 0), 0)
  expect_equal(mean(eztt[G[[1]]$genes, "First3.p"] <= 0.05 & eztt[G[[1]]$genes, "First3.logFC"] > 0), 0.1)
  
  #same behavior from fry & mroast
  expect_equal(rcn.f[,grep("^Prop", colnames(rcn.f))], rcn.m[,grep("^Prop", colnames(rcn.m))])
})

test_that("one pwy", {
  rcf1 <- roast_contrasts(object=M, G=G[1], feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="fry")
  expect_equal(rownames(rcf1), "pwy1")
  rcm1 <- roast_contrasts(object=M, G=G[1], feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="mroast")
  expect_equal(rownames(rcm1), "pwy1")
})

test_that("clean names & write", {
  unlink("bad_names_fry", recursive = TRUE) #in case it already exists
  
  G[[1]]$name <-"^path.way_1"
  G[[2]]$name <- "LPT2"
  expect_error(rcn.f <- roast_contrasts(object=M, G=G, feat.tab=eztt, grp=grp, contrast.v = contr.v, fun="fry", 
                                        name = "bad_names"), NA)
  unlink("bad_names_fry", recursive = TRUE)
})
  