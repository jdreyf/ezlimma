context("ezmediate")

# currently gives error without covariates
# takes a few sec -- worth it.
test_that("barfield", {
  prop.sig.size <- ezlimma:::sim_barfield(med.fcn = ezmediate, b1t2.v=0, nsim = 10, sims=30)
  expect_lte(prop.sig.size[1, 1], 0.2)
  
  prop.sig.size2 <- ezlimma:::sim_barfield(med.fcn = ezmediate, b1t2.v=0, nsim = 10, sims=30)
  expect_lte(prop.sig.size[1, 1], 0.2)
  
  prop.sig.power <- ezlimma:::sim_barfield(med.fcn = ezmediate, b1t2.v=0.39, nsim = 10, sims=30)
  expect_gte(prop.sig.power[1, 1], 0.2)
})
