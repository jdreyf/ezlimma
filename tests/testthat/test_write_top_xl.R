context("write top xl")

test_that("returned df", {
  fl <- lapply(G, FUN=function(x) x$genes)
  names(fl) <- lapply(G, FUN=function(x) x$name)
  wlx <- write_top_xl(pwy.tab=rcn.f, feat.lst=fl, feat.tab=eztt)
  
  expect_equal(grep("=HYPERLINK(", wlx[,1], fixed = TRUE), 1:nrow(wlx))
  expect_equal(wlx[,-1], rcn.f)
})