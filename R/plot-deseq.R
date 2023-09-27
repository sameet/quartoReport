#' Get the plotting data frame.
#'
#' @param res_df Results data frame for one comparison
#' @param thresh Threshold to determine significance. Default 0.05 or less.
#' @param label_n label these many genes.
#'
#' @return l a list of data frames that can be used to make plots.
# #' @export
#'
#' @examples
#' \dontrun{
#' use_l <- make_plot_df(res_df)
#' }
make_plot_df <- function(res_df, thresh = 0.05, label_n = 30){
  comparison <- extract_comparisons(res_df)

  res_df |>
    dplyr::filter(!is.na(padj)) |>
    dplyr::mutate(significant = ifelse(padj <= thresh, "yes", "no")) |>
    dplyr::mutate(updown = ifelse(log2FoldChange > 0, "up", "down")) |>
    dplyr::mutate(updown = factor(updown, levels = c("up", "down"))) -> plot_df

  op_l <- list(
    plot_df = plot_df,
    name    = comparison,
    label_n = label_n
  )

  op_l
}

#' Make a volcano plot from the result data-frame.
#'
#' @param res_df Results data frame.
#' @param ... Other parameters for the make plot df function.
#'
#' @return volcano_p a volcano plot, a ggplot2 graph/object
#' @export
#'
#' @examples
#' \dontrun{
#' volcano_p <- make_volcano(res_df)
#' }
make_volcano <- function(res_df, ...) {
  use_l <- make_plot_df(res_df, ...)

  use_l$plot_df |>
    dplyr::filter(significant == "yes") |>
    dplyr::arrange(desc(abs(log2FoldChange))) |>
    dplyr::slice_head(n = use_l$label_n) -> volcano_annot_df

  use_l$plot_df |>
    ggplot2::ggplot(ggplot2::aes(x = log2FoldChange,
                                 y = -10*log10(padj))) +
    ggplot2::geom_point(size = 0.2, col = "grey90", alpha = 0.1) -> p

  p + ggplot2::geom_point(data = use_l$plot_df |>
                            dplyr::filter(significant == "yes"),
                          ggplot2::aes(col = updown)) +
    ggrepel::geom_text_repel(data = volcano_annot_df,
                             ggplot2::aes(label = gene_id)) -> p
  p +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = use_l$name,
                  x = "log2(Fold Change)",
                  color = "Up/Down Expression",
                  y = "-10 x log10(padj)") +
    ggplot2::scale_color_manual(values = c("#2211FF", "#FF1122")) -> volcano_p

  volcano_p
}


#' My theme for plotting
#'
#' @return themed plot
#' @export
#'
#' @examples
#' \dontrun{
#' p + my_rnaseq_theme()
#' }
my_rnaseq_theme <- function(){
  ggplot2::theme(
    panel.border = ggplot2::element_rect(color = "#0C0C0C", fill = NA, linetype = 4, linewidth = 0.2),
    panel.background = ggplot2::element_rect(fill = "#F6F6F9"),
    axis.text = ggplot2::element_text(color = "#0B0B0B"),
    legend.background = ggplot2::element_rect(fill = "#F6F6F9", linetype = 3, linewidth = 0.25)
  )
}

#' Title
#'
#' @param res_df Results data frame for one comparison
#' @param vs_data Variance stabilized count data for the experiment
#' @param meta_df The meta-data frame
#' @param thresh Threshold to determne significance. Default padj <= 0.05
#' @param label_n Number of genes to plot.
#'
#' @return l a list of information to make plots
# #' @export
#'
#' @examples
#' \dontrun{
#' use_l <- make_box_plot_df(res_df, vs_data, meta_df)
#' }
make_box_plot_df <- function(res_df, vs_data, meta_df, thresh = 0.05, label_n = 30) {
  use_l <- make_plot_df(res_df, thresh = thresh, label_n = label_n)

  conditions <- unlist(strsplit(use_l$name, " -- "))

  use_samples <- meta_df |>
    dplyr::filter(condition %in% conditions) |>
    dplyr::select(sample, condition)

  use_l$plot_df |>
    dplyr::filter(significant == "yes") |>
    dplyr::arrange(desc(abs(log2FoldChange))) |>
    dplyr::slice_head(n = label_n) |>
    dplyr::arrange(log2FoldChange) |>
    dplyr::pull(gene_id) -> genes

  vs_data@assays@data[[1]] |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "gene_id") |>
    dplyr::filter(gene_id %in% genes) |>
    tidyr::pivot_longer(cols = !c("gene_id"), names_to = "sample", values_to = "expression") |>
    dplyr::inner_join(use_samples) |>
    dplyr::mutate(gene_id = factor(gene_id, levels = genes)) -> box_p_df

  op_l <- list(
    box_p_df = box_p_df,
    name = use_l$name
  )

  op_l
}

#' Make a single panel boxplot
#'
#' @param res_df Results data frame for one comparison
#' @param ... Other parameters to get the boxplot data-frame
#'
#' @return boxplot a ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' box_p <- make_box_single(res_df, vs_data, meta_df)
#' }
make_box_single <- function(res_df, ...) {
  use_l <- make_box_plot_df(res_df, ...)

  use_l$box_p_df |>
    ggplot2::ggplot(ggplot2::aes(x = gene_id, y = expression, fill = condition)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge2()) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   text = ggplot2::element_text(color = "#0F0F0F")) +
    ggplot2::labs(title = use_l$name,
         x = "Gene ID",
         y = "Normalized Expression",
         fill = "Condition(s)") +
    ggplot2::scale_fill_manual(values = c("#2222FF", "#FF2222")) +
    my_rnaseq_theme() -> p_box_single
  p_box_single
}


#' Create multipanel boxplot one gene per panel
#'
#' @param res_df Single result data frame
#' @param ... Other parameters to make_box_plot_df function
#'
#' @return multi-faceted-plot a multifacted boxplot ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' make_box_multi <- make_box_multi(res_df)
#' }
make_box_multi <- function(res_df, ...) {
  use_l <- make_box_plot_df(res_df, ...)

  use_l$box_p_df |>
    ggplot2::ggplot(ggplot2::aes(x = gene_id, y = expression, fill = condition)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge2(width = 0.95)) +
    ggplot2::facet_wrap(~ gene_id, scale = "free_x") +
    ggplot2::scale_fill_manual(values = c("#2222FF", "#FF2222")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 7),
                   strip.background = ggplot2::element_rect(fill = "#CBCBCB"),
                   axis.text.x.bottom = ggplot2::element_blank()) +
    ggplot2::labs(title = use_l$name,
         x = "",
         y = "Normalized Expression",
         fill = "Condition(s)") +
    my_rnaseq_theme() -> p_box_multi

  p_box_multi
}

#' Extract compairson from gven result file.
#'
#' @param res_df result data frame
#'
#' @return str comparison
#'
extract_comparisons <- function(res_df) {
  res_df |>
    dplyr::pull(comparison) |>
    unique() |>
    as.character()
}

#' Make df for an upSet plot
#'
#' @param dds A DESeqExperiment object
#' @param meta_df A meta-data data frame.
#' @param contrasts_df A contrasts data frame if specific contrasts in mind.  Can be null in which case all possible combinations will be used for plotting the graph.
#' @param thresh a threshold for padj value (default 0.05)
#'
#' @return upset_df A data frame that can be used to make an upSet plot
#' @export
#'
#' @examples
#' \dontrun{
#' upset_df <- make_upset_df(dds, meta_df, contrasts_df)
#' }
make_upset_df <- function(dds, meta_df, contrasts_df = NULL, thresh = 0.05) {
  all_res_df <- get_results_from_dds(dds = dds, df = contrasts_df)
  all_res_df |>
    dplyr::filter(padj <= 0.05) |>
    dplyr::mutate(comp_updown = ifelse(log2FoldChange > 0,
                                       paste(comparison, "up", sep = "_"),
                                       paste(comparison, "down", sep = "_"))) |>
    dplyr::group_by(gene_id) |>
    dplyr::summarize(comp_updown = list(comp_updown)) -> upset_df

  upset_df
}

#' Create an upSet plot
#'
#' @param upset_df A data frame with necessary list column. Generated by make_upset_df function call.
#'
#' @return upset_plot A ggplot2 object with upset parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' upset_plot <- make_upset_plot(upset_df) # when you have created the upset_df
#' upset_plot <- make_upset_plot(make_upset_df(dds, meta_df, contrasts_df, thresh))
#' # when you want to get to it directly.
#' }
make_upset_plot <- function(upset_df) {
  upset_df |>
    ggplot2::ggplot(ggplot2::aes(x = comp_updown)) +
    ggplot2::geom_bar() +
    ggplot2::geom_text(stat = "count",
                       ggplot2::aes(label = ggplot2::after_stat(count)),
                       vjust = 1.1, color = "#F0F0F0", size = 2) +
    ggplot2::labs(title = "Combined Gene Visualization",
                  x = "Comparison and Up/Down Status",
                  y = "Frequency (No. of Significant Genes)") +
    ggplot2::theme_minimal() +
    my_rnaseq_theme() +
    ggupset::scale_x_upset(order_by = "degree") +
    ggupset::theme_combmatrix(
      combmatrix.panel.point.color.empty = "grey90",
      combmatrix.panel.striped_background.color.one = "darkgrey",
      combmatrix.panel.striped_background.color.two = "white"
    ) -> upset_plot
  upset_plot
}

#' Title
#'
#' @param res_df Result data frame for a single comparison
#' @param meta_df Meta-data used in the analysis
#' @param vs_data Normalized gene expression after vst transformation
#' @param thresh Threshold to call significantly differentially expressed genes.
#'
#' @return op_l A list of objects that can be further manipulated.
#' @export
#'
#' @examples
#' \dontrun{
#' op_l <- make_heatmap(res_df, meta_df, vs_data) # will save a heatmap.png in the current directory
#' }
make_heatmap <- function(res_df, meta_df, vs_data, thresh = 0.05) {
  # Make a heatmap for a single comparison.  Include all the significant genes.

  # get the significant genes list from the results
  res_df |>
    dplyr::filter(padj <= thresh) |>
    dplyr::pull(gene_id) -> sig_genes

  # get the samples belonging to the conditions
  conds <- extract_comparisons(res_df = res_df) |>
    strsplit(" -- ") |>
    unlist()
  samples <- meta_df |>
    dplyr::filter(condition %in% conds)

  ofn <- paste(paste0(conds, collapse = "--"), "_sig_hm.png", sep = "")

  # get normalized expression counts for significant genes.
  vs_data@assays@data[[1]] |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "gene_id") |>
    dplyr::filter(gene_id %in% sig_genes) |>
    tidyr::pivot_longer(cols = !c("gene_id"), names_to = "sample", values_to = "expression") |>
    dplyr::inner_join(samples) |>
    tidyr::pivot_wider(id_cols = gene_id, names_from = sample, values_from = expression) -> hm_df

  # get annotation color dr
  col_annot_df <- meta_df |>
    dplyr::filter(condition %in% conds) |>
    dplyr::select(condition)

  my_colors <- list(
    condition = c("#F02211", "#1122F0")
  )
  names(my_colors$condition) <- conds

  # make the heatmap
  hm_df %>%
    tibble::column_to_rownames(var = "gene_id") %>%
    pheatmap::pheatmap(color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdGy")))(100),
                       border_color = NA,
                       scale = "row",
                       cluster_rows = TRUE,
                       cluster_cols = TRUE, main = "Heatmap (Significant Genes)",
                       annotation_col = col_annot_df,
                       annotation_colors = my_colors,
                       show_rownames = ifelse(nrow(.) <= 30, TRUE, FALSE),
                       silent = T) -> hm
  op_l <- list(
    hm_df = hm_df,
    col_annot_df = col_annot_df,
    hm = hm,
    ofn = ofn
  )

  op_l
}

#' save function for heatmap
#'
#' @param use_l Output from the make heatmap function call
#' @param use_dir Directory where you want to heatmap saved.
#'
#' @return returns nothing write the heatmap to the given folder with appropriate name.
#' @export
#'
#' @examples
#' \dontrun{
#' save_hm_pdf(use_l, "~/Desktop")
#' }
save_hm_pdf <- function(use_l, use_dir) {
  pdf(file = file.path(use_dir, gsub("png", "pdf", use_l$ofn)))
  grid::grid.newpage()
  grid::grid.draw(use_l$hm$gtable)
  dev.off()
}

#' save function for heatmap (png)
#'
#' @param use_l Output from the make heatmap function call
#' @param use_dir Directory where you want to heatmap saved.
#'
#' @return returns nothing write the heatmap to the given folder with appropriate name.
#' @export
#'
#' @examples
#' \dontrun{
#' save_hm_pdf(use_l, "~/Desktop")
#' }
save_hm_png <- function(use_l, use_dir) {
  png(file = file.path(use_dir, use_l$ofn))
  grid::grid.newpage()
  grid::grid.draw(use_l$hm$gtable)
  dev.off()
}
