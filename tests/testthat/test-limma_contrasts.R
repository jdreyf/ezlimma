context("limma_contrasts")

test_that("one contrast", {
  # p-values are sorted
  eztt1 <- limma_contrasts(M, grp = grp, contrast.v = contr.v[1])
  expect_equal(eztt1$First3.p, sort(eztt1$First3.p))
  
  contr.tmp <- c(First3="First3Arrays")
  eztt2 <- limma_contrasts(object=M, contrast.v = contr.tmp, design=design)
  expect_equal(eztt2$First3.p, sort(eztt2$First3.p))
})

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

test_that("error messages for duplicate or empty rownames", {
  M2 <- M
  rownames(M2)[2]  <- "gene1"
  expect_error(limma_contrasts(M2, grp = grp, contrast.v = contr.v))
  
  M3 <- M
  rownames(M3)[2] <- ""
  expect_error(limma_contrasts(M3, grp = grp, contrast.v = contr.v))
})

test_that("only one feature", {
  expect_error(limma_contrasts(object=M[1,], grp = grp, contrast.v = contr.v, add.means = FALSE))
  expect_silent(res <- limma_contrasts(object=M[1,,drop=FALSE], grp = grp, contrast.v = contr.v, add.means = FALSE))
  expect_equal(res[1, "Last3vsFirst3.logFC"], eztt["gene1", "Last3vsFirst3.logFC"])
  expect_silent(res <- limma_contrasts(object=M[1,,drop=FALSE], grp = grp, contrast.v = contr.v))
  expect_equal(res[1, "Last3vsFirst3.logFC"], eztt["gene1", "Last3vsFirst3.logFC"])
  expect_error(res <- limma_contrasts(object=data.matrix(t(M[1,])), grp = grp, contrast.v = contr.v))
  expect_equal(res[1, "Last3vsFirst3.logFC"], eztt["gene1", "Last3vsFirst3.logFC"])
})

test_that("!moderated", {
  expect_error(limma_contrasts(M, grp = grp, contrast.v = contr.v, trend=TRUE, moderated = FALSE))
  expect_error(limma_contrasts(M, grp = grp, contrast.v = contr.v, treat.lfc = 2, moderated = FALSE))
  
  eztt2 <- limma_contrasts(object=M, grp = grp, contrast.v = contr.v, moderated = FALSE,
                           cols=c("logFC", "t", "P.Value", "adj.P.Val"))
  g1.t <- t.test(x=M["gene1", grp=="Last3"], y=M["gene1", grp=="First3"], var.equal = TRUE)
  expect_equal(setNames(eztt2["gene1", "Last3vsFirst3.t"], nm="t"), g1.t$statistic)
  expect_equal(eztt2["gene1", "Last3vsFirst3.p"], g1.t$p.value)
  expect_equal(eztt2["gene1", "Last3vsFirst3.logFC"], 
               unname(g1.t$estimate["mean of x"]-g1.t$estimate["mean of y"]))
})

test_that("ndups & spacing", {
  expect_error(limma_contrasts(M, grp = grp, contrast.v = contr.v, block=rep(1:3, times=2), ndups=2))
  # need cor if ndups=2 or !is.null(block) 
  expect_error(limma_contrasts(M, grp = grp, contrast.v = contr.v, ndups=2))
  
  eztt.nodups <- limma_contrasts(object=M, grp = grp, contrast.v = contr.v, ndups = 1, spacing=1)
  expect_true(all.equal(eztt, eztt.nodups))
  
  eztt.dups <- limma_contrasts(object=M, grp = grp, contrast.v = contr.v, ndups = 2, correlation = 0.3)
  eztt.dups2 <- limma_contrasts(object=M, grp = grp, contrast.v = contr.v, ndups = 2, correlation = 0.6)
  expect_equal(nrow(eztt.dups), 50)
  expect_equal(dim(eztt.dups), dim(eztt.dups2))
  expect_true(!all(eztt.dups == eztt.dups2))
  # gene numbers in eztt.dups are all odds from 1:99:2
  expect_equal(sort(as.numeric(sub("gene", "", rownames(eztt.dups)))), seq(1, 99, by=2))
  expect_equal(eztt.dups["gene1", "First3.avg"], mean(M[1:2, 1:3]))
  expect_equal(eztt.dups["gene1", "Last3.avg"], mean(M[1:2, 4:6]))
  expect_true("gene1" %in% rownames(eztt.dups)[1:10])
  expect_lt(eztt.dups["gene1", "Last3vsFirst3.logFC"], 0)
  expect_gt(eztt.dups["gene1", "First3.avg"], 0.5)
  expect_lt(abs(mean(eztt.dups[, "Last3.avg"])), 0.1)
})