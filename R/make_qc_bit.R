#' Create the part to plot the QC file
#'
#' @param fn name of the rnaseq metrics file if available.
#'
#' @return str A string that will correctly plot the QC of the rnaseq metrics.
#' @export
#'
#' @examples
#' \dontrun{
#' qc_bit <- make_qc_bit()
#' }
make_qc_bit <- function(fn = params$metrics) {
  if(is.null(fn)) {
    qc_bit <- stringr::str_glue("
## QC Report
We do not have a metrics file available so we will do not have a QC section.
                                ")
  } else {
    qc_bit <- stringr::str_glue("
## QC Report
Usually (for most standard model organisms) a metrics file is generated.
This file helpfully named as `rnaseq-metrics.txt` is part of the output folder.
The metrics file reports on a number of metrics about the alignment for each sample.
Here we are showing a handful of those parameters for a quick assessment on each of the samples.
These metrics are generally created using the `PICARD` suite of programs.

```{{r}}
#| label: tab-showQC
#| tab-cap: QC parameters and their values/interpretation.

metrics_df <- read_delim(params$metrics, \"\\t\", show_col_types = FALSE)
metrics_df %>%
  dplyr::mutate(`% Unique` = gsub(\"%\", \"\", `% Unique`)) %>%
  dplyr::mutate(`% Unique` = as.numeric(`% Unique`)) -> metrics_df

metrics_df %>%
  dplyr::select(1, 3, 4, 5, 6, 7, 9, 10, 11, 12) %>%
  knitr::kable(booktabs = TRUE, linesep = \"\") %>%
  kableExtra::kable_paper(\"hover\", full_width = TRUE) %>%
  kableExtra::column_spec(1, color = ifelse(metrics_df$`% Unique` > 85, \"green\", \"red\"))

```
                                ")
  }
  qc_bit
}
