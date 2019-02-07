context("write top xl")

test_that("returned df", {
  wtx <- ezlimma:::write_top_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt)
  
  expect_equal(grep("=HYPERLINK(", wtx[,1], fixed = TRUE), 1:nrow(wtx))
  expect_equal(wtx[,-1], rcn.f)
  expect_equal(as.character(wtx[1,1]), '=HYPERLINK("pathways/pwy1.csv","pwy1")')
})

test_that("write out file", {
  ezlimma:::write_top_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt, name="text_wtx")
  expect_true(file.exists("text_wtx/text_wtx.xlsx"))
  unlink("text_wtx", recursive = TRUE)
})