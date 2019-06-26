context("write_gmt")

# also tested in test_read_gmt

test_that("write bad gmt.lst", {
  expect_error(write_gmt(G[[1]]))
})