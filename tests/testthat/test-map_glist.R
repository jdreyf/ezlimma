context("map glist")

test_that("helper G", {
  ann <- data.frame(sym=paste0("gene", 1:30, " /// gene", 30:1), stringsAsFactors = FALSE)
  rownames(ann) <- paste0("e", 1:30)
  mg <- map_glist(G=G, annot = ann, sep.str = " /// ", symbol.col = "sym")
  expect_true(all(mg[[1]]$genes == c(paste0("e", 1:10), paste0("e", 21:30))))
  
  mg2 <- map_glist(G=G, annot = ann, sep.str = " /// ", symbol.col = 1)
  expect_true(all(mg[[1]]$genes == mg2[[1]]$genes))
})