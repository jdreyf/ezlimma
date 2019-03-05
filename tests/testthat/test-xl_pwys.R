context("top xl")

test_that("returned df", {
  rf <- rcn.f
  rownames(rf)[2] <- paste0(rownames(rf)[2], ".")
  tx1 <- xl_pwys(pwy.tab=rf)
  expect_equal(grep("=HYPERLINK(", tx1[,1], fixed = TRUE), 1:nrow(tx1))
  expect_equal(rownames(tx1)[2], "pwy2_")
  expect_equal(as.character(tx1[2,1]), '=HYPERLINK("pathways/pwy2_.csv","pwy2.")')
})
