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

test_that("clean filenames", {
  expect_equal(clean_filenames("a"), "a")
  expect_equal(clean_filenames("."), "_")
  nm.test <- c("~", "`", "sulfate/sulfite met", "hi**", "hellohello", "hello world.", "CON", "coM2", "con.", "CON", "CON")
  res <- clean_filenames(nm.test)
  expect_false(any(duplicated(res)))
  expect_equal(res, c('_', '__', 'sulfate_sulfite_met', 'hi__', 'hellohello', "hello_world_",
                               'CON_', 'coM2_', 'con_', 'CON__', 'CON___'))
})