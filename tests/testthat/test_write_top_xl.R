context("write top xl")

test_that("returned df & written out file", {

  unlink("test_wtx", recursive = TRUE) #in case it already exists
  
  wtx <- write_top_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt, name="test_wtx")
  expect_equal(grep("=HYPERLINK(", wtx[,1], fixed = TRUE), 1:nrow(wtx))
  expect_equal(wtx[,-1], rcn.f)
  expect_equal(as.character(wtx[1,1]), '=HYPERLINK("pathways/pwy1.csv","pwy1")')
  unlink("test_wtx", recursive = TRUE)
  
  wtx1 <- write_top_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt, name="test_wtx", n.toptabs = 1)
  expect_equal(grep("=HYPERLINK(", wtx1[,1], fixed = TRUE), 1)
  unlink("test_wtx", recursive = TRUE)

  write_top_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt, name="test_wtx", n.toptabs = 1)
  expect_true(file.exists("test_wtx/test_wtx.xlsx"))
  expect_true(file.exists("test_wtx/pathways/pwy1.csv"))
  expect_false(file.exists("test_wtx/pathways/pwy2.csv"))
  unlink("test_wtx", recursive = TRUE)
})