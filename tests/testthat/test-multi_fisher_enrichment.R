context("multi fisher enrich")

test_that("non-sim multi fisher", {
  eztt.o <- eztt[order(eztt$Last3vsFirst3.p),]
  sig.set <- list(top=rownames(eztt.o)[1:50])
  fe.lst <- fisher_enrichment(sig.set = sig.set, G=G, feat.tab = eztt, return.lst = TRUE)
  fe <- fe.lst[[1]][order(fe.lst[[1]]$top.p),]
  
  ss2 <- list(bottom=rownames(eztt.o)[51:100])
  fe2 <- fisher_enrichment(sig.set = ss2, G=G, feat.tab = eztt)
  
  mfe <- multi_fisher_enrichment(sig.sets = c(sig.set, ss2), G=G, feat.tab = eztt)
  
  expect_equal(length(mfe), 2)
  expect_equal(length(mfe[[2]]), 3)
  expect_equal(mfe[[1]]["pwy2", "top.p"], fe["pwy2", "top.p"])
  expect_equal(mfe[[1]]["pwy1", "bottom.p"], fe2["pwy1", "bottom.p"])
  for (pwy.tmp in colnames(fe.lst$gene.membership)){
    mfe.tmp <- mfe$gene.membership[[pwy.tmp]][, "top", drop=FALSE]
    colnames(mfe.tmp) <- pwy.tmp
    expect_equal(fe.lst$gene.membership[rownames(mfe.tmp), pwy.tmp, drop=FALSE], mfe.tmp)
  }
})

