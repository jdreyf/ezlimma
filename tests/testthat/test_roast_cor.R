context("roast cor")

test_that("roast_cor weights", {
  #one less weight since have an NA
  set.seed(0)
  rcr.mw <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="mroast", weights=1:6)
  expect_equal(mean(rcr.mw$p==rcr.m$p), 0)
  #colnames of res don't start with .
  expect_equal(length(grep("^\\.", colnames(rcr.mw))), 0)
  rcr.fw <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="fry", weights=1:6)
  expect_equal(mean(rcr.fw$p==rcr.f$p), 0)
  
  #gene precision weights have effect
  rcr.fw <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="fry", weights=1:nrow(M))
  expect_equal(mean(rcr.f$p==rcr.fw$p), 0)
  set.seed(0)
  rcr.mw <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="mroast", weights=1:nrow(M))
  expect_equal(mean(rcr.m$p==rcr.mw$p), 0)
  
  #weights in EList
  rcr.fe <- roast_cor(el, G=G, stats.tab=eztt, pheno=pheno.v, fun="fry")
  expect_equal(mean(rcr.f$p==rcr.fe$p), 0)
  set.seed(0)
  rcr.me <- roast_cor(el, G=G, stats.tab=eztt, pheno=pheno.v, fun="mroast")
  expect_equal(mean(rcr.f$p==rcr.me$p), 0)
  
  #suppress object$weights
  expect_warning(rcr.fw <- roast_cor(el, G=G, stats.tab=eztt, pheno=pheno.v, fun="fry", weights=el$weights))
  expect_equal(mean(rcr.fe$p==rcr.fw$p), 1)
  set.seed(0)
  expect_warning(rcr.mw <- roast_cor(el, G=G, stats.tab=eztt, pheno=pheno.v, fun="mroast", weights=el$weights))
  expect_equal(mean(rcr.me$p==rcr.mw$p), 1)
  
  #pheno has NA
  expect_warning(rcr.mw2 <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno2, fun="mroast", weights=1:6))
  expect_equal(mean(rcr.mw$p==rcr.mw2$p), 0)
})

test_that("roast_cor one sided testing", {
  tmp <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="fry")
  tmp2 <- roast_cor(object=M, G=G, stats.tab=eztt, pheno=pheno.v, fun="fry", alternative = "less")
  expect_equal(tmp["pwy1", "p"]/2, tmp2["pwy1", "p"])
  expect_equal(tmp["pwy2", "p"]/2, tmp2["pwy2", "p"])
  expect_equal(tmp["pwy3", "p"]/2, tmp2["pwy3", "p"])
  #no mixed columns
  expect_equal(length(grep("Mixed", colnames(tmp2))), 0)
  
  set.seed(0)
  tmp3 <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="mroast")
  set.seed(0)
  tmp4 <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="mroast", alt="less")
  expect_lt(abs(tmp3["pwy1", "p"]/2 - tmp4["pwy1", "p"]), 0.001)
  expect_lt(abs(tmp3["pwy2", "p"]/2 - tmp4["pwy2", "p"]), 0.001)
  expect_lt(abs(tmp3["pwy3", "p"]/2 - tmp4["pwy3", "p"]), 0.001)
})

test_that("roast_cor trend has effect", {
  rc.res3 <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="fry", weights=1:6)
  rc.res3t <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="fry", weights=1:6, trend = TRUE)
  p.cols <- grep("\\.p", colnames(rc.res3))
  expect_equal(mean(rc.res3[,p.cols] == rc.res3t[,p.cols]), 0)
  
  set.seed(0)
  rc.res3 <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="mroast", weights=1:6)
  set.seed(0)
  rc.res3t <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="mroast", weights=1:6, trend = TRUE)
  p.cols <- grep("\\.p", colnames(rc.res3))
  expect_equal(mean(rc.res3[,p.cols] == rc.res3t[,p.cols]), 0)
})

test_that("p adjust", {
  expect_warning(roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="fry", adjust.method = "BY"))
  mr.by <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="mroast", adjust.method = "BY")
  expect_equal(colnames(mr.by)[ncol(mr.by)], "Mixed.BY")
})