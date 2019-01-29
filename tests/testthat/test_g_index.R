context('g_index')

G2 <- G
G2[[4]] <- list(name="pwy4", description=NA, genes=paste0("gene", 33:50))
G2[[5]] <- list(name="pwy5", description=NA, genes=paste0("gene", 31:32))

out.2 <- list(pwy1 = paste0('gene',1:10), pwy2 = paste0('gene',11:20), pwy3 = paste0('gene',21:30))

test_that('filters out the right gene sets', {
  expect_equal(g_index(G = G2, object = M, min.nfeats = 3, max.nfeats = 10), out.2)
  expect_error(g_index(G = G2, object = M, min.nfeats = 50, max.nfeats = 100))
})