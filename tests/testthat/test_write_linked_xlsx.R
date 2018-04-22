context("write_linked_xlsx")

tmp <- tempfile()
tmp.dir <- sub("\\\\file.+", "", tmp)
owd <- getwd()

setup({
  setwd(tmp.dir)
  library(xlsx)
})

test_that("file gets written", {
  rcn.f <- roast_contrasts(M, G=G, stats.tab=eztt, grp=grp, contrasts.v = contr.v, fun="fry", name="test")
  fry.file <- paste0(tmp.dir, "\\test_fry\\test_fry.xlsx")
  expect_true(file.exists(fry.file))
})

#can't detach xlsx
teardown({
  fry.dir <- paste0(tmp, "\\test_fry")
  unlink(fry.dir)
  setwd(owd)
})