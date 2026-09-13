context("top xl")

test_that("returned df", {
  rf <- rcn.f
  rownames(rf)[2] <- paste0(rownames(rf)[2], ".")
  tx1 <- xl_pwys(pwy.tab=rf)
  fmls <- xl_formula_text(tx1[,1])
  expect_equal(grep("=HYPERLINK(", fmls, fixed = TRUE), 1:nrow(tx1))
  expect_equal(rownames(tx1)[2], "pwy2_")
  expect_equal(fmls[2], '=HYPERLINK("pathways/pwy2_.csv","pwy2.")')
})
