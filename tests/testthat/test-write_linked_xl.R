context("write linked xl")

test_that("returned df & written out file", {
  wlx <- write_linked_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt, name="test_wlx")
  expect_equal(grep("=HYPERLINK(", wlx[,1], fixed = TRUE), 1:nrow(wlx))
  expect_equal(wlx[,-1], rcn.f)
  expect_equal(as.character(wlx[1,1]), '=HYPERLINK("pathways/pwy1.csv","pwy1")')
  expect_true(file.exists("test_wlx/test_wlx.xlsx"))
  expect_true(file.exists("test_wlx/pathways/pwy1.csv"))
  
  expect_silent(wlxr <- read_linked_xl("test_wlx/test_wlx.xlsx"))
  expect_equal(rcn.f, wlxr)
  
  rf <- rcn.f
  rownames(rf)[2] <- paste0(rownames(rf)[2], ".")
  rownames(rf)[3] <- "REACTOME_INHIBITION_OF_THE_PROTEOLYTIC_ACTIVITY_OF_APC_C_REQUIRED_FOR_THE_ONSET_OF_ANAPHASE_BY_MITOTIC_SPINDLE_CHECKPOINT_COMPONENTS"
  fl2 <- fl
  names(fl2)[2:3] <- rownames(rf)[2:3]
  wlx1 <- write_linked_xl(pwy.tab=rf, feat.lst=fl2, feat.tab=eztt, name="test_wlx", pwy.nchar = 100)
  expect_equal(grep("_", wlx1[,1]), 2:3)
  expect_true(file.exists("test_wlx/pathways/pwy2_.csv"))
  expect_true(file.exists("test_wlx/pathways/REACTOME_INHIBITION_OF_THE_PROTEOLYTIC_ACTIVITY_OF_APC_C_REQUIRED_FOR_THE_ONSET_OF_ANAPHASE_BY_MITOT.csv"))
  
  expect_silent(wlxr <- read_linked_xl("test_wlx/test_wlx.xlsx"))
  expect_equal(rf, wlxr)
})

teardown({
  unlink("test_wlx", recursive = TRUE)
})