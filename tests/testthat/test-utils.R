context("utils")

test_that("df_signif doesn't corrupt inputs", {
  test.df <- data.frame("num1-2"=1:4, num2=5:8, char1=letters[1:4], char2=letters[5:8], check.names=FALSE)
  expect_equal(df_signif(test.df[,2:4]), test.df[,2:4])
  expect_equal(df_signif(test.df[,1:3]), test.df[,1:3])
  expect_equal(df_signif(test.df[,3:4]), test.df[,3:4])
  expect_equal(df_signif(test.df[,1:2]), test.df[,1:2])
})

test_that("logfc2fc", {
  expect_equal(logfc2fc(2), 4)
  expect_equal(logfc2fc(0), 1)
  expect_equal(logfc2fc(-2), -4)
})

test_that("combine_pvals", {
  tab <- matrix(cbind(p=(1:10)/10, fdr=(1:10)/10), ncol=2, dimnames=list(letters[1:10], c("p", "fdr")))
  expect_equal(names(combine_pvalues(tab)), rownames(tab))
  tab2 <- data.frame(foo.p=(1:9)/9, bar.p=(9:1)/9)
  cp <- as.numeric(c('1', '0.626422003782904', '0.593346207732144', '0.581518654612405', '0.578313236039981', 
    '0.581518654612405', '0.593346207732144', '0.626422003782904', '1'))
  expect_equal(combine_pvalues(tab2), cp)
  #ensure grep for pv.cols not confused by "pheno.cor"
  tab2 <- data.frame(pheno.cor=(-2:2)/3, pheno.p=(1:5)/5)
  expect_equal(combine_pvalues(tab2), tab2[,2])
  #NULL = 0 rows
  expect_error(combine_pvalues(NULL))
})

test_that("clean filenames", {
  expect_equal(clean_filenames("a"), "a")
  expect_equal(clean_filenames("."), "_")
  nm.test <- c("~", "`", "sulfate/sulfite met", "hi**", "hellohello", "hello world.", "CON", "coM2", "con.", "CON", "CON")
  res <- clean_filenames(nm.test)
  expect_false(any(duplicated(res)))
  expect_equal(res, c('_', '__', 'sulfate_sulfite_met', 'hi__', 'hellohello', "hello_world_",
                               'CON_', 'coM2_', 'con_', 'CON__', 'CON___'))
})