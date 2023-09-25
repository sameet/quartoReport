#' Get the plotting data frame.
#'
#' @param res_df Results data frame for one comparison
#' @param thresh Threshold to determine significance. Default 0.05 or less.
#' @param label_n label these many genes.
#'
#' @return l a list of data frames that can be used to make plots.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_df <- make_plot_df(res_df)
#' }
make_plot_df <- function(res_df, thresh = 0.05, label_n = 30){
  comparison <- extract_comparison(res_df)

  res_df |>
    dplyr::filter(!is.na(padj)) |>
    dplyr::mutate(significant = ifelse(padj <= thresh, "yes", "no")) |>
    dplyr::mutate(updown = ifelse(log2FoldChange > 0, "up", "down")) -> plot_df

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
#'
#' @return volcano_p a volcano plot, a ggplot2 graph/object
#' @export
#'
#' @examples
#' \dontrun{
#' volcano_p <- make_volcano(res_df)
#' }
make_volcano <- function(res_df, ...) {
  use_l = make_plot_df(res_df, ...)

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
         y = "-10 x log10(padj)") -> volcano_p

  volcano_p
}

#' Extract compairson from gven result file.
#'
#' @param res_df
#'
#' @return str comparison
#'
extract_comparison <- function(res_df) {
  res_df |>
    dplyr::pull(comparison) |>
    unique() |>
    as.character()
}
