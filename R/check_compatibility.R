#' Check compatibility within the data
#'
#' @param count_mat A count matrix
#' @param meta_df A data frame generated from the sample-sheet.
#'
#' @return compatibility a logical value TRUE for compatible FALSE otherwise.
#' @export
#'
#' @examples
#' \dontrun{
#' count_mat <- read_count_mat("gene-count-matrix.csv")
#' meta_df <- read_meta_data("sample-sheet.txt")
#' check_compatibility(count_mat, meta_df)
#' }
check_compatibility <- function(count_mat, meta_df) {
  # check if the rownames of the count_mat match perfectly with meta_df
  if(all(colnames(count_mat) == rownames(meta_df))) TRUE
  else FALSE
}
