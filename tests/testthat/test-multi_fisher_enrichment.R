context("multi fisher enrich")

test_that("non-sim multi fisher", {
  eztt.o <- eztt[order(eztt$Last3vsFirst3.p),]
  sig.set.tmp <- list(top=rownames(eztt.o)[1:50])
  
  g2 <- G
  g2[[4]] <- list(name="pwy4", description=NA, genes=c("gene24", "gene25", "gene26"))
  
  fe <- fisher_enrichment(sig.set = sig.set.tmp, G=g2, feat.tab = eztt)
  fe <- fe[order(fe$top.p),]
  
  ss2 <- list(bottom=rownames(eztt.o)[51:100])
  fe2 <- fisher_enrichment(sig.set = ss2, G=g2, feat.tab = eztt)
  
  mfe <- multi_fisher_enrichment(sig.sets = c(sig.set.tmp, ss2), G=g2, feat.tab = eztt[, 1:3])
  
  expect_equal(nrow(mfe[[1]]), 4)
  expect_equal(mfe[[1]]["pwy4", "NGenes"], 3)
  expect_equal(ncol(mfe[[1]]), 7)
  expect_equal(mfe[[1]]["pwy2", "top.p"], fe["pwy2", "top.p"])
  expect_equal(mfe[[1]]["pwy1", "bottom.p"], fe2["pwy1", "bottom.p"])
  
  expect_equal(eztt.o[, "First3.p"], mfe[[2]][rownames(eztt.o), "First3.p"])
  expect_equal(mfe[[2]][, "top"], as.numeric(rownames(mfe[[2]]) %in% sig.set.tmp[[1]]))
  
  expect_equal(sum(is.na(rownames(mfe$pwy.stats))), 0)
  expect_equal(sum(is.na(rownames(mfe$feat.tab))), 0)
})