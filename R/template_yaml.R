#' Create opening YAML for a report
#'
#' @param title Title for the report
#' @param author Author for the report
#' @param email Email of the author for the report
#'
#' @return opening_yaml string to be used to construct the report qmd
#' @export
#'
#' @examples
#' \dontrun{
#' opnening_yaml <- template_yaml()
#' }
template_yaml <- function(title = "RNA-Seq Report",
                          author = "Sameet",
                          email = "sameet.mehta@yale.edu") {
  opening_yaml <- stringr::str_glue("
---
title: \"{title}\"
author: \"{author}\"
email: \"{email}\"
date: \"`r Sys.Date()`\"
format:
  pdf:
    coloredlinks: true
  html:
    code-fold: true
    code-summary: \"Show Code\"
  ipynb:
    markdown-headings: atx
    prefer-html: true
execute:
  warning: false
  echo: true
params:
  count_mat: \"gene_count_matrix.csv\"
  sample_sheet: \"sample-sheet.txt\"
  contrasts: null
  metrics: null
  outputs: \"report_output_dir\"
---
  ")
  opening_yaml
}
