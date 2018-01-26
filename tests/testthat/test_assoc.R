library("ezlimma")
context("assoc")

#example from limma::contrasts.fit
set.seed(42)
M <- matrix(rnorm(100*6, sd=0.3), nrow=100, ncol=6)
dimnames(M) <- list(paste0("gene", 1:nrow(M)), paste0("sample", 1:ncol(M)))
design <- cbind(First3Arrays=c(1,1,1,0,0,0), Last3Arrays=c(0,0,0,1,1,1))
grp <- rep(c("First3", "Last3"), each=3)
M[1, 1:3] <- M[1, 1:3] + 2

pheno.v <- rnorm(ncol(M))
pheno.v[1:3] <- pheno.v[1:3]-1
pheno2 <- pheno.v
pheno2[1] <- NA

contr.v <- c(First3="First3", Last3="Last3", Last3vsFirst3="Last3-First3")
eztt <- limma_contrasts(M, grp = grp, contrasts.v = contr.v)
eztt <- eztt[order(eztt$Last3vsFirst3.p),]

G <- list(list(name="pwy1", description=NA, genes=paste0("gene", 1:10)),
          list(name="pwy2", description=NA, genes=paste0("gene", 11:20)),
          list(name="pwy3", description=NA, genes=paste0("gene", 21:30)))

#######################################################################################################################
##gene-level
#######################################################################################################################
test_that("limma_contrasts", {
  fit <- lmFit(M, design=design)
  #  Would like to consider original two estimates plus difference between first 3 and last 3 arrays
  contrast.matrix <- cbind(First3=c(1,0),Last3=c(0,1),Last3vsFirst3=c(-1,1))
  fit2 <- contrasts.fit(fit,contrast.matrix)
  fit2 <- eBayes(fit2)
  toptab <- topTable(fit2, num=Inf, coef=3)
  
  eztt.ss <- eztt[,c("Last3vsFirst3.logFC", "Last3vsFirst3.p", "Last3vsFirst3.FDR")]
  col.nms <- c("logFC", "P.Value", "adj.P.Val")
  colnames(eztt.ss) <- col.nms
  expect_equal(eztt.ss, toptab[,col.nms])
  
  #give only p-val column
  eztt.p <- limma_contrasts(M, grp = grp, contrasts.v = contr.v, cols = "P.Value", add.means = FALSE)
  expect_equal(eztt.p, eztt[rownames(eztt.p), grep("\\.p$", colnames(eztt))])
})

test_that("contr_names", {
  contr.v <- c("First3")
  eztt <- limma_contrasts(M, grp = grp, contrasts.v = contr.v)
  expect_equal(length(grep("^\\.$", colnames(eztt))), 0)
})

test_that("ezcor", {
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

test_that("limma_cor", {
  design <- model.matrix(~1+pheno.v)
  fit <- lmFit(M, design)
  fit2 <- eBayes(fit)
  toptab <- topTable(fit2, coef=2, num=Inf)
  res2 <- limma_cor(M, phenotype = pheno.v)
  res2 <- res2[rownames(toptab),]
  # expect_equal(rownames(res2), rownames(toptab))
  expect_equal(res2$p, toptab$P.Value)
})

test_that("multi_cor", {
  res.ez <- ezcor(M, pheno2, method="spearman", reorder.rows = FALSE)
  res.mc <- multi_cor(M, cbind(a=pheno.v, b=pheno2), method="spearman", reorder.rows = FALSE)
  expect_equal(res.ez[,"p"], res.mc[,"b.p"])
  
  res.lm <- data.matrix(limma_cor(M, pheno2, reorder.rows = FALSE))
  res.mc <- multi_cor(M, cbind(a=pheno.v, b=pheno2), method="limma", reorder.rows = FALSE)
  expect_equal(res.lm[,"p"], res.mc[,"b.p"])
})

test_that("lmFit_weights", {
  wts <- (1:ncol(M))/ncol(M)
  fit.w <- lmFit(M, design=design, weights=wts)
  fit.aw <- lmFit(M, design=design, array.weights=wts)
  fit0 <- lmFit(M, design=design)
  expect_equal(fit0, fit.aw)
  #identical(fit0, fit.aw)
})

#######################################################################################################################
##pwy-level
#######################################################################################################################
test_that("roast_contrasts", {
  rc.res1 <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry", weights=1:6)
  expect_equal(rownames(rc.res1)[1], "pwy1")
  rc.res2 <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast", weights=1:6)
  expect_equal(rownames(rc.res2)[1], "pwy1")
})

test_that("roast_cor", {
  #one less weight since have an NA
  rc.res3 <- roast_cor(M, G=G, stats.tab=eztt, pheno=pheno.v, fun="mroast", weights=1:6)
  expect_equal(rownames(rc.res3)[1], "pwy1")
  #colnames of res don't start with .
  expect_equal(length(grep("^\\.", colnames(rc.res3))), 0)
})

test_that("roast_contrasts one sided testing", {
  tmp <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry")
  tmp2 <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry", alternative = "less")
  expect_equal(1-tmp["pwy1", "First3.p"]/2, tmp2["pwy1", "First3.p"])
  expect_equal(tmp["pwy2", "First3.p"]/2, tmp2["pwy2", "First3.p"])
  #no mixed columns
  expect_equal(length(grep("Mixed", colnames(tmp2))), 0)
  
  tmp3 <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast")
  tmp4 <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="mroast", alt="less")
  expect_lt(abs(1-tmp3["pwy1", "First3.p"]/2 - tmp4["pwy1", "First3.p"]), 0.02)
  expect_lt(abs(tmp3["pwy2", "First3.p"]/2 - tmp4["pwy2", "First3.p"]), 0.02)
  expect_lt(abs(tmp3["pwy3", "First3.p"]/2 - tmp4["pwy3", "First3.p"]), 0.02)
})