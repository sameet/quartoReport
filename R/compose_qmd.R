#' Compose the Quarto markdown file.
#'
#' @param ofn Name of the output file.  A temporary file will be generated
#'
#' @return ofn Path to the output file.
#' @export
#'
#' @examples
#' \dontrun{
#' ofn <- compose_qmd()
#' }
compose_qmd <- function(ofn = NULL) {
  # compose the complete qmd file.
  if(is.null(ofn)) {
    ofn <- tempfile(pattern = "rnaseq-report-", fileext = ".qmd")
  }
  message(paste("Using ", ofn, sep = ""))

  # gather all components
  opening_yaml <- template_yaml(author = "Sameet Mehta", email = "sameet.mehta@yale.edu")
  setup_part <- make_setup_qmd()
  qc_bit <- make_qc_bit(fn = "params$metrics")
  analysis1 <- make_analysis()

  sink(ofn)
  print(opening_yaml)
  print(setup_part)
  print(qc_bit)
  print(analysis1)
  sink()

  ofn
}
