context("fc2ratio")

test_that("fc2ratio", {
  expect_equal(fc2ratio(2), (2^8 - 2^7)/2^7)
  expect_equal(fc2ratio(3), 2)
  expect_equal(fc2ratio(4), (2^8 - 2^6)/2^6)
  expect_equal(fc2ratio(4.5), 3.5)

  #use example of logFC with log2 expr ~= 6
  expect_equal(fc2ratio(-1.5), (2^(6 - log2(1.5)) - 2^6)/2^6)
  expect_equal(fc2ratio(FC = -2), (2^5 - 2^6)/2^6)
  expect_equal(fc2ratio(-8), (2^4 - 2^7)/2^7)
  expect_equal(fc2ratio(-4), (2^4 - 2^6)/2^6)
})
