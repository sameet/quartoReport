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
  meta_fn: \"sample-sheet.txt\"
  contrasts_fn: null
  metrics: null
  outputs: \"report_output_dir\"
  use_threshold: 0.05
---
  ")
  opening_yaml
}

#' Make setup part of the qmd file
#'
#' @return setup_qmd string that holds some text to setup the analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' setup_qmd <- make_setup_qmd()
#' }
make_setup_qmd <- function() {
  setup_qmd <- stringr::str_glue("

## Introduction

This is analysis report for differential gene expression analysis.
The data is described in the later sections.
In this section we are setting up the requirements for the analysis to happen correcty.
We will import packages and functions.

```{{r}}
#| label: \"setup\"

library(quartoReport)
# use_thresh <- 0.05 # global threshold to use as alpha in the analysis
```

Please note that this package is currently under active development.
This status may change in the future.
Currently this package is **not** ready for production use.
USE AT YOUR OWN RISK.

                                 ")
  # message(class(setup_qmd))
  setup_qmd
}
