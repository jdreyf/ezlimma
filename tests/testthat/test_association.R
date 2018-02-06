library("ezlimma")
context("association")

#example from limma::contrasts.fit
set.seed(0)
M <- matrix(rnorm(100*6, sd=0.3), nrow=100, ncol=6)
dimnames(M) <- list(paste0("gene", 1:nrow(M)), paste0("sample", 1:ncol(M)))
design <- cbind(First3Arrays=c(1,1,1,0,0,0), Last3Arrays=c(0,0,0,1,1,1))
grp <- rep(c("First3", "Last3"), each=3)
M[1, 1:3] <- M[1, 1:3] + 2

el <- new("EList")
el$E <- M
ww <- matrix(rexp(n=nrow(M)*ncol(M)), ncol=ncol(M), nrow=nrow(M))
el$weights <- ww

pheno.v <- rnorm(ncol(M))
pheno.v[1:3] <- pheno.v[1:3]-1
pheno2 <- pheno.v
pheno2[1] <- NA

contr.v <- c(First3="First3", Last3="Last3", Last3vsFirst3="Last3-First3")
eztt <- limma_contrasts(M, grp = grp, contrasts.v = contr.v)
lc <- limma_cor(M, phenotype = pheno.v)

G <- list(list(name="pwy1", description=NA, genes=paste0("gene", 1:10)),
          list(name="pwy2", description=NA, genes=paste0("gene", 11:20)),
          list(name="pwy3", description=NA, genes=paste0("gene", 21:30)))

rcn.f <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry")
set.seed(0)
rcn.m <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast")
rcr.f <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="fry")
set.seed(0)
rcr.m <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="mroast")

#######################################################################################################################
##gene-level
#######################################################################################################################
test_that("limma_contrasts matches topTable(eBayes(contrasts.fit(lmfit(M))))", {
  fit <- lmFit(M, design=design)
  #  Would like to consider original two estimates plus difference between first 3 and last 3 arrays
  contrast.matrix <- cbind(First3=c(1,0),Last3=c(0,1),Last3vsFirst3=c(-1,1))
  fit2 <- contrasts.fit(fit,contrast.matrix)
  fit2 <- eBayes(fit2)
  toptab <- topTable(fit2, num=Inf, coef=3)
  
  eztt.ss <- eztt[,c("Last3vsFirst3.logFC", "Last3vsFirst3.p", "Last3vsFirst3.FDR")]
  col.nms <- c("logFC", "P.Value", "adj.P.Val")
  colnames(eztt.ss) <- col.nms
  expect_equal(eztt.ss, toptab[rownames(eztt.ss), col.nms])
  
  #give only p-val column
  eztt.p <- limma_contrasts(M, grp = grp, contrasts.v = contr.v, cols = "P.Value", add.means = FALSE)
  expect_equal(eztt.p, eztt[rownames(eztt.p), grep("\\.p$", colnames(eztt))])
})

test_that("limma_contrasts trend has effect", {
  #give only p-val column
  ezttt.p <- limma_contrasts(M, grp = grp, contrasts.v = contr.v, cols = "P.Value", add.means = FALSE, trend=TRUE)
  #most p-values should be different
  expect_equal(mean(ezttt.p == eztt[rownames(ezttt.p), grep("\\.p$", colnames(eztt))]), 0)
})

test_that("empty contr_names don't create cols that start with '.'", {
  contr.v <- c("First3")
  eztt.nonm <- limma_contrasts(M, grp = grp, contrasts.v = contr.v)
  expect_equal(length(grep("^\\.$", colnames(eztt.nonm))), 0)
})

test_that("ezcor with different methods matches cor", {
  #make separate apply calls to ensure no bugs
  m.r <- apply(M, 1, function(v){
    cor(v, pheno.v)
  })
  m.p <- apply(M, 1, function(v){
    cor.test(v, pheno.v)$p.value
  })
  res1 <- ezcor(M, pheno.v, method="pearson", reorder.rows = FALSE)
  expect_equal(res1[,"cor"], m.r)
  expect_equal(res1[,"p"], m.p)
  
  #make separate apply calls to ensure no bugs
  m.r <- apply(M, 1, function(v){
    cor(v, pheno.v, method = "kendall")
  })
  m.p <- apply(M, 1, function(v){
    cor.test(v, pheno.v, method = "kendall")$p.value
  })
  res1 <- ezcor(M, pheno.v, method="kendall", reorder.rows = FALSE)
  expect_equal(res1[,"tau"], m.r)
  expect_equal(res1[,"p"], m.p)
})

test_that("limma_cor matches topTable(eBayes(lmFit(M)))", {
  design <- model.matrix(~1+pheno.v)
  fit <- lmFit(M, design)
  fit2 <- eBayes(fit)
  toptab <- topTable(fit2, coef=2, num=Inf)
  # expect_equal(rownames(res2), rownames(toptab))
  expect_equal(lc[rownames(toptab), "p"], toptab$P.Value)
})

test_that("limma_cor trend has effect", {
  lc2 <- limma_cor(M, phenotype = pheno.v, trend=TRUE)
  expect_equal(mean(lc2$p == lc[rownames(lc2), "p"]), 0)
})

test_that("multi_cor matches ezcor & limma_cor", {
  res.ez <- ezcor(M, pheno2, method="spearman", reorder.rows = FALSE)
  res.mc <- multi_cor(M, cbind(a=pheno.v, b=pheno2), method="spearman", reorder.rows = FALSE)
  expect_equal(res.ez[,"p"], res.mc[,"b.p"])
  
  res.lm <- data.matrix(limma_cor(M, pheno2, reorder.rows = FALSE))
  res.mc <- multi_cor(M, cbind(a=pheno.v, b=pheno2), method="limma", reorder.rows = FALSE)
  expect_equal(res.lm[,"p"], res.mc[,"b.p"])
})

test_that("lmFit array weights only affect if given in 'weights'", {
  wts <- (1:ncol(M))/ncol(M)
  fit.aw <- lmFit(M, design=design, array.weights=wts)
  fit0 <- lmFit(M, design=design)
  expect_equal(fit0, fit.aw)
  
  fit.w <- lmFit(M, design=design, weights=wts)
  expect_false(identical(fit0, fit.w))
  
  fit.gw <- lmFit(M, design=design, weights=1:nrow(M))
  expect_false(identical(fit0, fit.w))
})

test_that("limma_contrasts weights", {
  eztt.w <- limma_contrasts(M, grp = grp, contrasts.v = contr.v, weights = 1:ncol(M))
  expect_equal(mean(eztt.w$First3.p==eztt$First3.p), 0)
  
  eztt.gw <- limma_contrasts(M, grp = grp, contrasts.v = contr.v, weights = 1:nrow(M))
  expect_equal(mean(eztt.gw$First3.p==eztt$First3.p), 0)
  
  #create EList object
  eztt.el <- limma_contrasts(el, grp = grp, contrasts.v = contr.v)
  expect_equal(mean(eztt.el$First3.p==eztt$First3.p), 0)
  
  #object$weights are being ignored
  expect_warning(eztt.elw <- limma_contrasts(el, grp = grp, contrasts.v = contr.v, weights = NULL))
  expect_equal(mean(eztt.elw$First3.p==eztt$First3.p), 1)
})

test_that("limma_cor weights", {
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
})

#######################################################################################################################
##pwy-level
#######################################################################################################################
test_that("roast_contrasts weights", {
  #aw have effect
  rcn.fw <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry", weights=1:6)
  expect_equal(mean(rcn.f$First3.p==rcn.fw$First3.p), 0)
  set.seed(0)
  rcn.mw <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast", weights=1:6)
  expect_equal(mean(rcn.m$First3.p==rcn.mw$First3.p), 0)
  
  #gene weights have effect
  expect_warning(rcn.fw <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry", gene.weights=1:nrow(M)))
  expect_equal(mean(rcn.f$First3.p==rcn.fw$First3.p), 1)
  set.seed(0)
  rcn.mw <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast", gene.weights=1:nrow(M))
  expect_equal(mean(rcn.m$First3.p==rcn.mw$First3.p), 0)
  
  rcn.fe <- roast_contrasts(el, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry")
  expect_equal(mean(rcn.f$First3.p==rcn.fe$First3.p), 0)
  set.seed(0)
  rcn.me <- roast_contrasts(el, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast")
  expect_equal(mean(rcn.m$First3.p==rcn.me$First3.p), 0)
  
  #suppress object$weights
  expect_warning(rcn.fw <- roast_contrasts(el, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry", weights=el$weights))
  expect_equal(mean(rcn.fw$First3.p==rcn.fe$First3.p), 1)
  set.seed(0)
  expect_warning(rcn.mw <- roast_contrasts(el, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast", weights=el$weights))
  expect_equal(mean(rcn.mw$First3.p==rcn.me$First3.p), 1)
})

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
})

test_that("roast_contrasts one sided testing", {
  tmp <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry")
  tmp2 <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry", alternative = "less")
  expect_equal(1-tmp["pwy1", "First3.p"]/2, tmp2["pwy1", "First3.p"])
  expect_equal(tmp["pwy2", "First3.p"]/2, tmp2["pwy2", "First3.p"])
  expect_equal(1-tmp["pwy3", "First3.p"]/2, tmp2["pwy3", "First3.p"])
  #no mixed columns
  expect_equal(length(grep("Mixed", colnames(tmp2))), 0)
  
  set.seed(42)
  tmp3 <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast")
  set.seed(42)
  tmp4 <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast", alt="less")
  expect_lt(abs(1-tmp3["pwy1", "First3.p"]/2 - tmp4["pwy1", "First3.p"]), 0.001)
  expect_lt(abs(tmp3["pwy2", "First3.p"]/2 - tmp4["pwy2", "First3.p"]), 0.001)
  expect_lt(abs(1-tmp3["pwy3", "First3.p"]/2 - tmp4["pwy3", "First3.p"]), 0.001)
})

test_that("roast_contrasts trend has effect", {
  rc.res1 <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry", weights=1:6)
  rc.res1t <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry", weights=1:6, trend = TRUE)
  expect_equal(mean(rc.res1$First3.p == rc.res1t$First3.p), 0)
  
  set.seed(0)
  rc.res2 <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast", weights=1:6)
  set.seed(0)
  rc.res2t <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast", weights=1:6, trend = TRUE)
  p.cols <- grep("\\.p", colnames(rc.res2))
  expect_equal(mean(rc.res2[,p.cols] == rc.res2t[,p.cols]), 0)
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