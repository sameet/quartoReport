test_that("read_count_mat(fn) without meta-data file stops", {
  expect_error(read_count_mat(testthat::test_path("..",
                                                    "test_data",
                                                    "gene_count_matrix_fixed.csv"),
                                meta_df = NULL))
})
test_that("read_count_mat() read a file and returns an integer matrix", {
  meta_df <- read_meta_data(testthat::test_path("..", "test_data", "sample-sheet.txt"))
  count_mat <- read_count_mat(testthat::test_path("..", "test_data", "gene_count_matrix_fixed.csv"),
                              meta_df)
  expect_gt(nrow(count_mat), 50000)
  # expect_(count_mat, "matrix")
})
