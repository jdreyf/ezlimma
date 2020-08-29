context("grep cols")

test_that("from helper", {
  expect_equal(grep_cols(tab=rcr.m, p.cols = "p")[1], "p")
  expect_equal(grep_cols(tab=rcr.m, p.cols = 5), "p")
  expect_error(grep_cols(tab=rcr.m[, 1:3], p.cols = "p"))
  expect_error(grep_cols(tab=rcr.m, stat.cols = "Down"))
  # ok to have one NA
  rcr.m[1, "p"] <- NA
  expect_silent(grep_cols(tab=rcr.m, p.cols = 5))
  # not ok to have full column be NA
  rcr.m[2:3, "p"] <- NA
  expect_error(grep_cols(tab=rcr.m, p.cols = 5))
  
  expect_equal(grep_cols(tab=lc, stat.cols = 4), "slope")
  expect_equal(grep_cols(tab=lc, stat.cols = "slope"), "slope")
})