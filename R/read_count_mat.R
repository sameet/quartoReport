#' Read in the counts
#'
#' @param fn path to a file name typically a csv file.
#' @param meta_df A meta-data frame.
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
read_count_mat <- function(fn, meta_df = NULL) {
  if(!fs::is_file(fn)) {
    stop("Need count file.")
  }

  if(is.null(meta_df)) {
    stop("Need a meta data data-frame")
  }

  if(!is.data.frame(meta_df)) {
    stop("Meta-data needs to be a data-frame")
  }

  fn |>
    readr::read_delim(",", show_col_types = FALSE) |>
    tidyr::pivot_longer(!tidyselect::starts_with("gene_id"),
                        names_to = "sample", values_to = "counts") |>
    dplyr::group_by(gene_id, sample) |>
    dplyr::summarize(fixed_counts = sum(counts,
                                        na.rm = TRUE)) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(id_cols = gene_id,
                        names_from = sample,
                        values_from = fixed_counts,
                        values_fill = 0) |>
    as.data.frame() |>
    tibble::column_to_rownames(var = "gene_id") |>
    as.matrix() -> count_mat
  count_mat <- count_mat[, rownames(meta_df)]
  count_mat
}
