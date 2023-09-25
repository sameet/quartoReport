#' Create DESeq2 object from counts and sample-sheet
#'
#' @param count_mat count matrix from the count data
#' @param meta_df meta-data data-frame from sample-sheet
#'
#' @return A DESeq2 object that can be used for further analysis
#' @export
#'
#' @examples
#' \dontrun{
#' dds <- make_deseq_object(count_mat, meta_df)
#' }
make_deseq_object <- function(count_mat, meta_df) {
  compatibility <- check_compatibility(count_mat = count_mat,
                                       meta_df = meta_df)
  if(!compatibility) stop("Data not compatible")

  dds <- suppressMessages(
    suppressWarnings(
      DESeq2::DESeqDataSetFromMatrix(countData = count_mat,
                                     colData = meta_df,
                                     design = ~condition,
                                     tidy = FALSE)
    )
  )

  dds <- DESeq2::DESeq(dds) # the actual calculation happen here
  dds
}

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

