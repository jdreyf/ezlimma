context("read_gmt")

tmp <- tempfile()
setup({
  g.tab <- rbind(G[[1]]$genes, G[[2]]$genes, G[[3]]$genes)
  g.tab <- cbind(c(G[[1]]$name, G[[2]]$name, G[[3]]$name), rep(NA, 3), g.tab)
  write.table(g.tab, tmp, col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE)
})

test_that("read_gmat", {
  rg <- read_gmt(tmp)
  expect_equal(names(rg), c(G[[1]]$name, G[[2]]$name, G[[3]]$name))
  expect_equal(rg[[2]]$genes, G[[2]]$genes)
})

teardown({
  unlink(tmp)
})