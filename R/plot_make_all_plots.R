#' Title
#'
#' @param res_df A result data frame for a single comparison
#' @param thresh Threshold to be used for calling significant genes
#' @param op_dir Output directory where the plots should be saved.
#' @param vs_data A vst stabilized data object.
#' @param meta_df A data frame containing sample metadata. Must have 'sample' and 'condition' columns.
#' @param label_n Number of genes to be lablled, or number of boxes in boxplot
#' @param interactive Logical, whether to convert ggplot objects to interactive plotly objects. Defaults to FALSE.
#'
#' @return plot_l A list of plots.
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate static plots
#' plot_l_static <- make_all_plots(res_df, vs_data, meta_df, thresh, op_dir)
#' # Generate interactive plots
#' plot_l_interactive <- make_all_plots(res_df, vs_data, meta_df, thresh, op_dir, interactive = TRUE)
#' # plot_l_interactive$volcano will now be a plotly object
#' # plot_l_interactive$bp_s will also be a plotly object
#' }
make_all_plots <- function(res_df, vs_data, meta_df, thresh = 0.05, label_n = 30, op_dir = NULL, interactive = FALSE) {

  if(is.null(op_dir)) {
    message("Graphical outputs will not be saved.")
  }

  # The following could be a function of its own, we may refactor it later.
  res_df |>
    dplyr::filter(padj <= thresh) |>
    dplyr::pull(gene_id) |>
    length() -> n_sig_genes

  if(n_sig_genes < label_n) {
    label_n <- n_sig_genes
  }

  vp <- make_volcano(res_df = res_df, thresh = thresh, label_n = label_n)
  bp_s <- make_box_single(res_df = res_df, vs_data = vs_data, meta_df = meta_df,
                          label_n = label_n, thresh = thresh)

  if (interactive) {
    vp <- to_interactive_plot(vp)
    bp_s <- to_interactive_plot(bp_s)
  }

  hm_l <- make_heatmap(res_df = res_df, meta_df = meta_df, vs_data = vs_data, thresh = thresh)
  plot_l <- list(
    volcano = vp,
    bp_s = bp_s,
    hm = hm_l$hm
  )

  res_df %>%
    extract_comparisons() %>%
    gsub(" -- ", "--", .) -> ofn_pref

  if(!is.null(op_dir)) {
    message(paste("Saving files to ", op_dir))

    ggplot2::ggsave(filename = file.path(op_dir, paste(ofn_pref, "-volcano.pdf", sep = "")),
                    plot = vp, width = 9, height = 9, units = "in")
    ggplot2::ggsave(filename = file.path(op_dir, paste(ofn_pref, "-boxplot_singlePanel.pdf", sep = "")),
                    plot = bp_s, width = 9, height = 9, units = "in")
    ggplot2::ggsave(filename = file.path(op_dir, paste(ofn_pref, "-sig_hm.pdf", sep = "")),
                    plot = ggplotify::as.ggplot(hm_l$hm))
  }

  plot_l
}

