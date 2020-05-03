context("fisher enrich")

# sim to test fisher_enrichment
test_that("sim_fisher", {
  alpha <- 0.25
  fet <- sim_fisher(G=G, grp=grp, feat.tab=eztt, nsim=100, alpha=alpha, verbose=FALSE)
  
  # sim also tests with stopifnot
  expect_lte(fet[1, 1], alpha)
  expect_gte(fet[1, 2], fet[1, 1])
})

test_that("non-sim_fisher", {
  eztt.o <- eztt[order(eztt$Last3vsFirst3.p),]
  sig.set <- list(top=rownames(eztt.o)[1:50])
  fe.lst <- fisher_enrichment(sig.set = sig.set, G=G, feat.tab = eztt, return.lst = TRUE)
  fe <- fe.lst[[1]][order(fe.lst[[1]]$top.p),]
  expect_equal(rownames(fe)[1], "pwy1")
  expect_equal(as.numeric(colSums(fe.lst[[2]])), fe.lst[[1]]$top.num)
  
  ss2 <- list(bottom=rownames(eztt.o)[51:100])
  fe2 <- fisher_enrichment(sig.set = ss2, G=G, feat.tab = eztt)
  expect_true("bottom.p" %in% colnames(fe2))
})
  