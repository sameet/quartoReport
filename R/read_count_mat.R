#' Read in the counts
#'
#' @param fn
#'
#' @return A matrix of counts of integers.
#' @export
#'
#' @examples
#' count_mat <- read_count_mat("genes_count_mat.csv")
read_count_mat <- function(fn) {
  if(!fs::is_file(fn)) {
    stop("Need count file.")
  }

  fn |>
    readr::read_delm(fn, ",") |>
    as.matrix()
}
