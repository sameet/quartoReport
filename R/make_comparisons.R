#' Create combinations of comparions from given condition vector
#'
#' @param v A character vector with at least 2 conditions
#' @param n Number of conditions to use in combination, default is 2
#'
#' @return combination_df A data frame of possible combinations
#' @export
#'
#' @examples
#' combination_df <- make_comparisons(c("Control", "Test1", "Test2"))
make_comparisons <- function(v, n = 2) {
  # expect a vector of strings, and we make a data frame of 2 conditions on each row.

  if(!is.character(v)) stop("The condition's vector should be a character")

  v |>
    utils::combn(n) |>
    t() |>
    as.data.frame() |>
    stats::setNames(c("c1", "c2")) -> comparison_df

  comparison_df
}
