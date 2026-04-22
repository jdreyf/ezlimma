test_that("works on tibbles", {
  v1 <- combine_pvalues(eztt)
  v2 <- tibble::as_tibble(eztt) |> combine_pvalues()
  v3 <- eztt |> tibble::rowid_to_column("gene") |> combine_pvalues()
  expect_equal(setNames(v1, nm=NULL), v2)
  expect_equal(v2, v3)
})
