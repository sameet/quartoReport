test_that("multiplication works", {
  expect_equal(class(read_meta_data(testthat::test_path("tests",
                                                      "test_data",
                                                      "sample-sheet.txt"))), "data.frame")
})
