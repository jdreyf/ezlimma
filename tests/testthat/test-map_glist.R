context("map glist")

test_that("helper G", {
  ann <- data.frame(sym=paste0("gene", 1:30, " /// gene", 30:1), stringsAsFactors = FALSE)
  rownames(ann) <- paste0("e", 1:30)
  mg <- map_glist(G=G, annot = ann, sep.str = " /// ", symbol.col = "sym")
  expect_true(all(mg[[1]]$genes == c(paste0("e", 1:10), paste0("e", 21:30))))
  
  mg2 <- map_glist(G=G, annot = ann, sep.str = " /// ", symbol.col = 1)
  expect_true(all(mg[[1]]$genes == mg2[[1]]$genes))
})

test_that("new G", {
  # want to convert "g1" to "gene1"
  G.tmp <- list(list(name="pwy1", description=NA, genes=paste0("g", 1:10)),
            list(name="pwy2", description=NA, genes=paste0("g", 11:20)),
            list(name="pwy3", description=NA, genes=paste0("g", 21:30)))
  
  annot <- data.frame(Gene=rownames(M), row.names=rownames(M))
  annot.tmp <- data.frame(annot, sym=paste0("g", 1:100))
  expect_equal(rownames(annot.tmp), rownames(M))
  annot.tmp$sym[1] <- "g1 __ g11"
  
  ecc <- map_glist(G=G.tmp, annot=annot.tmp, symbol.col = "sym", sep.str = " __ ")
  expect_equal(length(ecc), 3)
  expect_equal(length(ecc[[1]][[3]]), 10)
  expect_equal(length(ecc[[2]][[3]]), 11)
  expect_equal(length(ecc[[3]][[3]]), 10)
  expect_true("gene1" %in% ecc[[1]][[3]])
  expect_true("gene1" %in% ecc[[2]][[3]])
  expect_false("gene1" %in% ecc[[3]][[3]])
  
  annot.tmp$sym[1] <- "g1|g11"
  ecb <- map_glist(G=G.tmp, annot=annot.tmp, symbol.col = "sym", sep.str = "|", fixed = TRUE)
  expect_true("gene1" %in% ecb[[1]][[3]])
  expect_true("gene1" %in% ecb[[2]][[3]])
  expect_false("gene1" %in% ecb[[3]][[3]])
})
