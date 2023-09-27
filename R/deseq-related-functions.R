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

#' Get results from a dds object
#'
#' @param dds An object of class DESeqExperiment
#' @param df A data frame with contrasts, one on each row. Each row contains only two columns, named c1, and c2
#'
#' @return Gives a data-frame for the results with those comparisons.
#' @export
#'
#' @examples
#' \dontrun{
#' df <- get_results_from_dds(dds, condition_df)
#' }
get_results_from_dds <- function(dds = dds,
                                 df = NULL) {
  if(is.null(df)){
    res_df <- get_all_results(dds)
  } else {
    lapply(seq_len(nrow(df)), function(i) {
      get_single_result(dds = dds, comp_df = df[1, ])
    }) %>%
      do.call("rbind", .) -> res_df
  }
  res_df
}


#' Get results for all contrasts.
#'
#' @param dds A DESeqExperiment object.
#'
#' @return all_res_df A data-frame with results for all the contrasts.
#' @export
#'
#' @examples
#' \dontrun{
#' all_res_df <- get_all_results(dds)
#' }
get_all_results <- function(dds){
  use_comp_df <- get_comp_from_dds(dds)
  lapply(seq_len(nrow(comp_df)), function(i) {
    df <- get_single_result(dds, use_comp_df[i, ])
  }) %>%
    do.call("rbind", .) -> all_res_df
  all_res_df
}

#' Extract possible comparisons from dds colData slot
#'
#' @param dds A DESeqDataSet Object
#'
#' @return comb_df a data frame with possible combination of comparisons
#'
#' @examples
#' \dontrun{
#' comb_df <- get_comp_from_dds(dds)
#' }
get_comp_from_dds <- function(dds) {
  SummarizedExperiment::colData(dds) |>
    as.data.frame() |>
    dplyr::pull(condition) |>
    as.character() |>
    unique() |>
    combn(2) |>
    t() |>
    as.data.frame() |>
    setNames(c("c1", "c2")) -> comp_df

  comp_df
}

#' Get single result from combination df, and dds object
#'
#' @param dds object of class DESeqDataSet
#' @param comp_df one row of a combination df
#'
#' @export
#' @return res_df A data frame with result comparing the two conditions.
#' @examples
#' \dontrun{
#' res_df <- get_single_result(dds, comp_df)
#' }
get_single_result <- function(dds, comp_df) {
  # the same function will allow returning results from specific contrasts file.
  cond_l <- comp_df |> as.list()
  cond_str <- paste0(unlist(cond_l), collapse = " -- ")
  res <- DESeq2::results(dds, contrast = c("condition", cond_l$c1, cond_l$c2))
  res <- res %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene_id") %>%
    dplyr::mutate(comparison = cond_str)

  res
}

#' Save result files.
#'
#' @param res_df Result data frame for a single result
#' @param op_dir Output directory
#' @param thresh Threshold to determine significance (default padj <= 0.05)
#'
#' @return No return value.
#' @export
#'
#' @examples
#' \dontrun{
#' save_results(res_df, op_dir = tempdir())
#' }
save_results <- function(res_df, op_dir, thresh = 0.05) {
  # save the result table and significant genes table for a comparison
  if(is.null(op_dir)) {
    stop("Needs an output directory to save the files")
  } else {
    message(paste("Saving files to :", op_dir))
  }

  res_df %>%
    extract_comparisons() %>%
    gsub(" -- ", "--", .) -> ofn_prefix

  n_sig_genes <- res_df |>
    dplyr::filter(padj <= thresh) |>
    dplyr::pull(gene_id) |>
    length()

  ofn_all <- file.path(op_dir, paste(ofn_prefix, "_all-results.txt", sep = ""))
  ofn_sig <-
  res_df |>
    readr::write_delim(ofn_all, delim = "\t")

  if(n_sig_genes > 0) {
    ofn_sig <- file.path(op_dir, paste(ofn_prefix, "_sig-results.txt", sep = ""))
    res_df |>
      dplyr::filter(padj <= thresh) |>
      readr::write_delim(ofn_sig, delim = "\t")
  } else {
    message("No significant gene detected for the comparison ", ofn_prefix)
  }
}
