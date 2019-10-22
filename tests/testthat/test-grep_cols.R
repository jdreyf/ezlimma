context("grep cols")

test_that("from helper", {
  expect_equal(grep_cols(tab=rcr.m, p.cols = "p")[1], "p")
  expect_equal(grep_cols(tab=rcr.m, p.cols = 5), "p")
  rcr.m[1, "p"] <- NA
  expect_error(grep_cols(tab=rcr.m, p.cols = 5))
  expect_error(grep_cols(tab=rcr.m, stat.cols = "Down"))
    
  expect_equal(grep_cols(tab=lc, stat.cols = 4), "slope")
  expect_equal(grep_cols(tab=lc, stat.cols = "slope"), "slope")
})