#' Read Meta data (sample sheet)
#'
#' @param fn Tab separated values file with samples and conditions
#'
#' @return meta_df a data.frame object.
#' @export
#'
#' @examples
#' \dontrun{
#' meta_df <- read_meta_data(file.path("tests/test_data/sample-sheet.txt"))
#' }
read_meta_data <- function(fn) {
  fn |>
    readr::read_delim("\t", show_col_types = FALSE) |>
    stats::setNames(c("sample", "condition")) |>
    as.data.frame() -> meta_df
  rownames(meta_df) <- meta_df$sample
  meta_df
}
