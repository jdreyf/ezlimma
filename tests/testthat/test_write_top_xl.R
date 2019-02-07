context("write top xl")

test_that("returned df", {
  wtx <- write_top_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt)
  
  expect_equal(grep("=HYPERLINK(", wtx[,1], fixed = TRUE), 1:nrow(wtx))
  expect_equal(wtx[,-1], rcn.f)
  expect_equal(as.character(wtx[1,1]), '=HYPERLINK("pathways/pwy1.csv","pwy1")')
  
  wtx1 <- write_top_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt, n.toptabs = 1)
  expect_equal(grep("=HYPERLINK(", wtx1[,1], fixed = TRUE), 1)
})

test_that("write out file", {
  unlink("text_wtx", recursive = TRUE) #in case it already exists
  
  write_top_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt, name="text_wtx", n.toptabs = 1)
  expect_true(file.exists("text_wtx/text_wtx.xlsx"))
  expect_true(file.exists("text_wtx/pathways/pwy1.csv"))
  expect_false(file.exists("text_wtx/pathways/pwy2.csv"))
  unlink("text_wtx", recursive = TRUE)
})