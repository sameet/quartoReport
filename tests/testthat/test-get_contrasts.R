test_that("get_contrasts(NULL, NULL) fails with an error", {
  expect_error(get_contrasts())
})

test_that("get_contrasts(fn = NULL, meta_df) works", {
  meta_df <- read_meta_data(testthat::test_path("..",
                                                "test_data",
                                                "sample-sheet.txt"))
  contrast_df <- get_contrasts(fn = NULL, meta_df)
  expect_equal(nrow(contrast_df), 6)
})
