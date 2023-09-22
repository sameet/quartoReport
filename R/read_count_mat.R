#' Read in the counts
#'
#' @param fn path to a file name typically a csv file.
#'
#' @return A matrix of counts of integers.
#' @export
#'
#' @examples
#' \dontrun{
#' library(here)
#' fn <- file.path("tests", "test_data", "gene_count_matrix_fixed.csv")
#' count_mat <- read_count_mat(fn)
#' }
read_count_mat <- function(fn) {
  if(!fs::is_file(fn)) {
    stop("Need count file.")
  }

  fn |>
    readr::read_delim(",") |>
    as.matrix()
}
