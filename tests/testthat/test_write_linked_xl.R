context("write linked xl")

test_that("returned df & written out file", {
  wlx <- write_linked_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt, name="test_wlx")
  expect_equal(grep("=HYPERLINK(", wlx[,1], fixed = TRUE), 1:nrow(wlx))
  expect_equal(wlx[,-1], rcn.f)
  expect_equal(as.character(wlx[1,1]), '=HYPERLINK("pathways/pwy1.csv","pwy1")')
  expect_true(file.exists("test_wlx/test_wlx.xlsx"))
  expect_true(file.exists("test_wlx/pathways/pwy1.csv"))
  
  rf <- rcn.f
  rownames(rf)[2] <- paste0(rownames(rf)[2], ".")
  fl2 <- fl
  names(fl2)[2] <- rownames(rf)[2]
  wlx1 <- write_linked_xl(pwy.tab=rf, feat.lst=fl2, feat.tab=eztt, name="test_wlx")
  expect_equal(grep("_", wlx1[,1]), 2)
  expect_true(file.exists("test_wlx/pathways/pwy2_.csv"))
})

teardown({
  unlink("test_wlx", recursive = TRUE)
})