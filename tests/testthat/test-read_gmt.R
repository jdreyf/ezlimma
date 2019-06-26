context("read_gmt")

tmp <- tempfile()

test_that("write & read gmt", {
  write_gmt(G, tmp)
  rg <- read_gmt(tmp)
  # only rg has names; G$description is NA, rg$description is "NA"
  expect_equal(names(rg), c(G[[1]]$name, G[[2]]$name, G[[3]]$name))
  expect_equal(rg[[2]]$genes, G[[2]]$genes)
})

teardown({
  unlink(tmp)
})