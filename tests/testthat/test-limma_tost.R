context("limma tost")

test_that("tost", {
  # gene 1
  lt <- limma_tost(object=M, grp = grp, contrast.v = contr.v[3], tost.lfc = log2(1.2))
  expect_gte(lt["gene1", "Last3vsFirst3_tost.p"], 0.9)
  expect_gte(lt["gene1", "Last3vsFirst3_tost.t"], max(lt[setdiff(rownames(lt), "gene1"), "Last3vsFirst3_tost.t"]))

  # only gene94
  # TOSTER::TOSTtwo.raw(m1=lt["gene94", "First3.avg"], m2=lt["gene94", "Last3.avg"], 
  #                     sd1=sd(M["gene94", grp=="First3"]), sd2=sd(M["gene94", grp=="Last3"]),
  #                     n1=3, n2=3, low_eqbound=-log2(1.2), high_eqbound=log2(1.2), alpha=0.05, var.equal = TRUE)
  # gives p = 0.0425!
  
  lt2 <- limma_tost(object=M["gene94",,drop=FALSE], grp = grp, contrast.v = contr.v[3], tost.lfc = log2(1.2), add.means = FALSE)
  expect_equal(as.numeric(signif(lt2[1, "Last3vsFirst3_tost.p"], 3)), 0.0425)
})