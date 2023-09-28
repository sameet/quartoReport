make_pca_obj <- function(df) {
  df %>%
    t() %>%
    stats::prcomp() %>%
    return()
}


get_sig_genes <- function(all_res_df, thresh = 0.05) {
  # get a list of singificant genes from all results data frame.
  all_res_df |>
    dplyr::filter(padj <= thresh) |>
    dplyr::pull(gene_id) |>
    unique()
}

#' Make a PCA plot from a DESeq2 object
#'
#' @param dds A DESeqExperiment object
#' @param vs_data A vst transformed data object
#' @param contrasts_df A data frame of contrasts of interest
#' @param thresh Threshold to call significant genes (default 0.05)
#'
#' @return pca_plot_l a list of ggplot2 objects.
#' @export
#'
#' @examples
#' \dontrun{
#' pca_plot_l <- make_pca_plot_l(dds, vs_data, contrasts_df)
#' }
make_pca_plot_l <- function(dds, vs_data, contrasts_df, thresh = 0.05) {
  # Create a PCA plot at given threshold.
  vs_data@assays@data[[1]] |>
    as.data.frame() -> all_vs_df

  dds |>
    SummarizedExperiment::colData() |>
    as.data.frame() |>
    dplyr::select(sample, condition) -> use_meta_df

  all_vs_df |>
    make_pca_obj() -> all_genes_pca

  all_res <- get_results_from_dds(dds = dds,
                                  df = contrasts_df)

  all_res |>
    get_sig_genes(thresh = thresh) -> sig_genes

  all_vs_df |>
    tibble::rownames_to_column(var = "gene_id") |>
    dplyr::filter(gene_id %in% sig_genes) |>
    tibble::column_to_rownames(var = "gene_id") |>
    make_pca_obj() -> sig_genes_pca

  p1 <- make_pca_plot(pca_obj = all_genes_pca, meta_df = use_meta_df, use_title = "PCA (all genes)")
  p2 <- make_pca_plot(pca_obj = sig_genes_pca, meta_df = use_meta_df, use_title = "PCA (sig-only)")
  p3 <- make_scree_plot(pca_obj = all_genes_pca)
  p4 <- make_scree_plot(pca_obj = sig_genes_pca)

  op_l <- list(
    all_pca = p1,
    sig_pca = p2,
    all_scree = p3,
    sig_scree = p4
  )
  op_l
}

#' Make PCA plot from a pca object
#'
#' @param pca_obj PCA Object
#' @param meta_df meta-data data frame. Can be a colData from dds
#' @param pc_1 First PC to be used in the scatter plot.
#' @param pc_2 Second PC to be used in the scatter plot.
#' @param use_title The title for the subplot
#'
#' @return p ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' p <- make_pca_plot(pca_obj, meta_df, use_title = "PCA Plot")
#' }
make_pca_plot <- function(pca_obj, meta_df, use_title, pc_1 = 1, pc_2 = 2) {
  # Make PCA scatter plot and color by the conditions of the samples.

  pc1var <- round(summary(pca_obj)$importance[2, pc_1] * 100, digits = 2)
  pc2var <- round(summary(pca_obj)$importance[2, pc_2] * 100, digits = 2)

  use_pc_1 = paste0("PC", pc_1, collapse = "")
  use_pc_2 = paste0("PC", pc_2, collapse = "")

  pca_obj$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample") %>%
    dplyr::inner_join(meta_df, by = "sample") %>%
    dplyr::mutate(condition = factor(condition, levels = unique(condition))) %>% #head()
    ggplot2::ggplot(ggplot2::aes(x = get(use_pc_1), y = get(use_pc_2),
                                 color = condition, label = sample)) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel() +
    ggplot2::labs(
      title = use_title,
      x = paste(pc1var, "% Captured ", use_pc_1, sep = ""),
      y = paste(pc2var, "% Captured ", use_pc_2, sep = "")
    ) +
    ggplot2::theme_minimal() +
    my_rnaseq_theme() -> p
  p
}

#' Scree plot to explain importance of the PCs calculated
#'
#' @param pca_obj An object generated with prcomp function call
#'
#' @return p ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' p <- make_scree_plot(pca_obj)
#' }
make_scree_plot <- function(pca_obj) {
  # Make a summary plot to show importance of the Principal Components

  summary(pca_obj)$importance %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "metric") %>%
    tidyr::pivot_longer(cols = !metric) %>%
    dplyr::mutate(name = factor(name, levels = paste("PC", 1:nrow(.), sep = ""))) %>%
    dplyr::filter(metric == "Proportion of Variance") %>%
    dplyr::mutate(value = 100*value) %>%
    ggplot2::ggplot(ggplot2::aes(x = name, y = value)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(
      x = "PC (name)",
      y = "[%] Contribution"
    ) +
    ggplot2::theme_bw() +
    my_rnaseq_theme() -> p
  p
}
