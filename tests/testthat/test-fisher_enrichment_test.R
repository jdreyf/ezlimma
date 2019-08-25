context("fisher")

test_that("sim", {
  alpha <- 0.25
  fet <- sim_fisher(G=G, grp=grp, feat.tab=eztt, nsim=100, alpha=alpha, verbose=FALSE)

  # sim also tests with stopifnot
  expect_lte(fet[1, 1], alpha)
  expect_gte(fet[1, 2], fet[1, 1])
})