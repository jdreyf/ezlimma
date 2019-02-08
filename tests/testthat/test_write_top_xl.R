context("write top xl")

test_that("returned df & written out file", {

  unlink("test_wtx", recursive = TRUE) #in case it already exists
  
  wtx <- write_top_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt, name="test_wtx")
  expect_equal(grep("=HYPERLINK(", wtx[,1], fixed = TRUE), 1:nrow(wtx))
  expect_equal(wtx[,-1], rcn.f)
  expect_equal(as.character(wtx[1,1]), '=HYPERLINK("pathways/pwy1.csv","pwy1")')
  expect_true(file.exists("test_wtx/test_wtx.xlsx"))
  expect_true(file.exists("test_wtx/pathways/pwy1.csv"))
  unlink("test_wtx", recursive = TRUE)
  
  rf <- rcn.f
  rownames(rf)[2] <- paste0(rownames(rf)[2], ".")
  fl2 <- fl
  names(fl2)[2] <- rownames(rf)[2]
  wtx1 <- write_top_xl(pwy.tab=rf, feat.lst=fl2, feat.tab=eztt, name="test_wtx")
  expect_equal(grep("_", wtx1[,1]), 2)
  unlink("test_wtx", recursive = TRUE)
})