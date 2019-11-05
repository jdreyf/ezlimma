context("clean filenames")

test_that("clean filenames", {
  expect_equal(clean_filenames("a"), "a")
  expect_equal(clean_filenames("."), "_")
  nm.test <- c("~", "`", "sulfate/sulfite met", "hi**", "hellohello", "hello world.", "CON", "coM2", "con.", "CON", "CON")
  res <- clean_filenames(nm.test)
  expect_false(any(duplicated(res)))
  expect_equal(res, c('_', '__', 'sulfate_sulfite_met', 'hi__', 'hellohello', "hello_world_",
                      'CON_', 'coM2_', 'con_', 'CON__', 'CON___'))
  
  # long nms
  nm <- "REACTOME_INHIBITION_OF_THE_PROTEOLYTIC_ACTIVITY_OF_APC_C_REQUIRED_FOR_THE_ONSET_OF_ANAPHASE_BY_MITOTIC_SPINDLE_CHECKPOINT_COMPONENTS"
  expect_equal(nchar(clean_filenames(nm, nm.nchar = 100)), 100)
})