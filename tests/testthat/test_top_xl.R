context("top xl")

test_that("returned df", {
  tx1 <- top_xl(pwy.tab=rcn.f, n.toptabs = 1)
  expect_equal(grep("=HYPERLINK(", tx1$xl[,1], fixed = TRUE), 1)
  expect_equal(tx1$pwy.csv.nms, "pwy1")
})
