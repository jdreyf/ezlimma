context("limma_contrasts")

test_that("lmFit array weights only affect if given in 'weights'", {
  wts <- (1:ncol(M))/ncol(M)
  fit.aw <- lmFit(M, design=design, array.weights=wts)
  fit0 <- lmFit(M, design=design)
  expect_equal(fit0, fit.aw)
  
  fit.w <- lmFit(M, design=design, weights=wts)
  expect_false(identical(fit0, fit.w))
  
  fit.gw <- lmFit(M, design=design, weights=1:nrow(M))
  expect_false(identical(fit0, fit.gw))
})

test_that("matches topTable(eBayes(contrasts.fit(lmfit(M))))", {
  fit <- lmFit(M, design=design)
  #  Would like to consider original two estimates plus difference between first 3 and last 3 arrays
  contrast.matrix <- cbind(First3=c(1,0),Last3=c(0,1),Last3vsFirst3=c(-1,1))
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  toptab <- topTable(fit2, num=Inf, coef=3)
  
  eztt.ss <- eztt[,c("Last3vsFirst3.logFC", "Last3vsFirst3.p", "Last3vsFirst3.FDR")]
  col.nms <- c("logFC", "P.Value", "adj.P.Val")
  colnames(eztt.ss) <- col.nms
  expect_equal(eztt.ss, toptab[rownames(eztt.ss), col.nms])
  
  #give only p-val column
  eztt.p <- limma_contrasts(M, grp = grp, contrast.v = contr.v, cols = "P.Value", add.means = FALSE)
  expect_equal(eztt.p, eztt[rownames(eztt.p), grep("\\.p$", colnames(eztt))])
})

test_that("trend has effect", {
  #give only p-val column
  ezttt.p <- limma_contrasts(M, grp = grp, contrast.v = contr.v, cols = "P.Value", add.means = FALSE, trend=TRUE)
  #most p-values should be different
  expect_equal(mean(ezttt.p == eztt[rownames(ezttt.p), grep("\\.p$", colnames(eztt))]), 0)
})

test_that("empty contr_names don't create cols that start with '.'", {
  contr.v <- c("First3")
  eztt.nonm <- limma_contrasts(M, grp = grp, contrast.v = contr.v)
  expect_equal(length(grep("^\\.$", colnames(eztt.nonm))), 0)
})

test_that("weights", {
  eztt.w <- limma_contrasts(M, grp = grp, contrast.v = contr.v, weights = 1:ncol(M))
  expect_equal(mean(eztt.w$First3.p==eztt$First3.p), 0)
  
  eztt.gw <- limma_contrasts(M, grp = grp, contrast.v = contr.v, weights = 1:nrow(M))
  expect_equal(mean(eztt.gw$First3.p==eztt$First3.p), 0)
  
  #create EList object
  eztt.el <- limma_contrasts(el, grp = grp, contrast.v = contr.v)
  expect_equal(mean(eztt.el$First3.p==eztt$First3.p), 0)
  
  #object$weights are being ignored
  expect_warning(eztt.elw <- limma_contrasts(el, grp = grp, contrast.v = contr.v, weights = NULL))
  expect_equal(mean(eztt.elw$First3.p==eztt$First3.p), 1)
})

test_that("with treat matches topTable(treat(contrasts.fit(lmfit(M))))", {
  fit <- lmFit(M, design=design)
  #  Would like to consider original two estimates plus difference between first 3 and last 3 arrays
  contrast.matrix <- c(Last3vsFirst3=c(1, -1))
  expect_warning(fit2 <- contrasts.fit(fit, contrasts=contrast.matrix))
  fit2 <- treat(fit2, lfc=log2(1.1))
  toptab <- topTable(fit2, num=Inf, coef=1, sort.by="p")
  
  #give only p-val column
  eztt.p <- limma_contrasts(M, grp = grp, contrast.v = contr.v[3], cols = "P.Value", add.means = FALSE, treat.lfc = log2(1.1))
  expect_equal(eztt.p[,1], toptab$P.Value)
})


test_that("error messages for duplicate or empty rownames",{
  M2 <- M
  rownames(M2)[2]  <- 'gene1'
  expect_error(limma_contrasts(M2, grp = grp, contrast.v = contr.v))
  
  M3 <- M
  rownames(M3)[2] <- ''
  expect_error(limma_contrasts(M3, grp = grp, contrast.v = contr.v))
})