#' Stabilize count data using vst
#'
#' @param dds a DESeq2 object
#'
#' @return vs_data a count object with vst transformed counts.
#' @export
#'
#' @examples
#' \dontrun{
#' vs_data <- stabilize_count_data(dds)
#' }
stabilize_count_data <- function(dds) {
  vs_data <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  vs_data
}
