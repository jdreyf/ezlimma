context("fisher")

test_that("sim_fisher", {
  alpha <- 0.25
  fet <- sim_fisher(G=G, grp=grp, feat.tab=eztt, nsim=100, alpha=alpha, verbose=FALSE)

  # sim also tests with stopifnot
  expect_lte(fet[1, 1], alpha)
  expect_gte(fet[1, 2], fet[1, 1])
})

test_that("non-sim_fisher", {
  eztt.o <- eztt[order(eztt$Last3vsFirst3.p),]
  sig.sets <- list(top=rownames(eztt.o)[1:50], bottom=rownames(eztt.o)[51:100])
  fe <- fisher_enrichment(sig.sets = sig.sets, G=G, feat.tab = eztt)
  fe <- fe[order(fe$top.p),]
  expect_true("bottom.p" %in% colnames(fe))
  expect_equal(rownames(fe)[1], "pwy1")
})