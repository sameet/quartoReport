test_that("read_meta_data(fn) yeilds a dataframe", {
  expect_equal(class(read_meta_data(testthat::test_path("..",
                                         "test_data",
                                         "sample-sheet.txt"))), "data.frame")
})
