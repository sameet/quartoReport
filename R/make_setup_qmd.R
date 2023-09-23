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
use_thresh <- 0.05 # global threshold to use as alpha in the analysis
```

Please note that this package is currently under active development.
This status may change in the future.
Currently this package is **not** ready for production use.
USE AT YOUR OWN RISK.

                                 ")
  # message(class(setup_qmd))
  setup_qmd
}
