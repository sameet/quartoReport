#' Get contrasts for analysis
#'
#' @param fn name of the contrasts file.
#' @param meta_df meta-data data-frame we use this if contrasts file is not provided.
#'
#' @return contrasts_df a contarsts data-frame that will be used to extract data from the DESeq2 object.
#' @export
#'
#' @examples
#' \dontrun{
#' contrasts_df <- get_contrasts("contrasts.txt", meta_df)
#' contarsts_df <- get_contrasts(meta_df = meta_df)
#' }
get_contrasts <- function(fn = NULL, meta_df = NULL) {
  if(is.null(fn)) {
    if(is.null(meta_df)) rlang::abort("Need either a contrasts file name or meta_df")

    else {
      message("We do not have a contrasts file.  We will use the meta-data to generate the contrasts")
      meta_df |>
        dplyr::pull(condition) |> # extracts a column and makes it into a vector
        unique() |>
        make_comparisons() -> contrasts_df
    }

  } else {

    fn |>
      readr::read_delim("\t", show_col_types = FALSE) |>
      stats::setNames(c("c1", "c2")) -> contrasts_df
  }
  contrasts_df
}
