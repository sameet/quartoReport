#' Get contrasts for analysis
#'
#' @param fn name of the contrasts file.
#' @param meta_df meta-data data-frame we use this if contrasts file is not provided.
#'
#' @return contrasts_df a contarsts data-frame that will be used to extract data from the DESeq2 object.
#' @export
#'
#' @examples
get_contrasts <- function(fn = NULL, meta_df) {
  if(fs::is_file(fn)) {
    fn |>
      readr::read_delim("\t", show_col_types = FALSE) |>
      stats::setNames(c("c1", "c2")) -> contrasts_df
  }
  if(is.null(fn)) {
    message("We do not have a contrasts file.  We will use the meta-data to generate the contrasts")
    meta_df |>
      pull(condition) |>
      unique() |>
      make_comparisons() -> contrasts_df
  }
  contrasts_df
}
