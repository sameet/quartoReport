#' Compose the Quarto markdown file.
#'
#' @param ofn Name of the output file.  A temporary file will be generated
#' @param contrasts_df The contrasts_df to determine number of contrasts in the data.
#' @param thresh The threshold for calculation of default 0.05
#' @param label_n Number of genes to label if possible default 30
#' @param metric_fn Name of the file with RNA-Seq alignment metrics.
#' @param ... Other parameters to pass to opening yaml function
#'
#' @return ofn Path to the output file.
#' @export
#'
#' @examples
#' \dontrun{
#' ofn <- compose_qmd()
#' }
compose_qmd <- function(ofn = NULL, contrasts_df, thresh = 0.05, label_n = 30, metric_fn = NULL, ...) {
  # compose the complete qmd file.
  if(is.null(ofn)) {
    ofn <- tempfile(pattern = "rnaseq-report", fileext = ".qmd")
  }
  message(paste("Using ", ofn, sep = ""))

  # gather all components
  opening_yaml <- template_yaml(author = "Sameet Mehta", email = "sameet.mehta@yale.edu", ...)
  setup_part <- make_setup_qmd()
  qc_bit <- make_qc_bit(fn = metric_fn)
  analysis1 <- make_analysis()
  overall_scope <- make_overall_scope_bit(contrasts_df)


  # lapply(seq_len(nrow(contrasts_df)), function(i) {
  #   my_ind <- i
  #   comparison_part <- make_single_comparison_bit(get_single_result(dds = dds, comp_df = contrasts_df[my_ind, ]),
  #                                                 thresh = 0.05, label_n = 30, comp_n = i)
  # }) %>%
  #   unlist() %>%
  #   paste0(collapse = "\n") -> comparison_part_all

  sink(ofn)
  print(opening_yaml)
  print(setup_part)
  print(qc_bit)
  print(analysis1)
  lapply(seq_len(nrow(contrasts_df)), function(i) {
    use_res_df <- get_single_result(dds, comp_df = contrasts_df[i, ])
    comparison_part <- make_single_comparison_bit(res_df = use_res_df,
                                                  thresh = 0.05, label_n = 30, comp_n = i)
    print(comparison_part)
  })
  print(overall_scope)
  sink()

  ofn
}


