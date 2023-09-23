test_that("make_comparisons(c('Control', 'Test1', 'Test2')) creates a data frame", {
  comparison_df <- make_comparisons(c("Control", "Test1", "Test2"))
  expect_equal(class(comparison_df), "data.frame")
  expect_equal(nrow(comparison_df), 3)
})

test_that("make_comparisons() Stops with an error", {
  expect_error(make_comparisons())
})
