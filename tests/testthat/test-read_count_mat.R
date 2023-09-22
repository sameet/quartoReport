test_that("read_count_mat() read a file and returns an integer matrix", {
  expect_gt(nrow(read_count_mat(testthat::test_path("..", "test_data", "gene_count_matrix_fixed.csv"))),
            50000)
})
