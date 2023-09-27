#' Save computed DESEq2 objects.
#'
#' @param obj Object to be saved, either a DESEq2 experiment or vst data
#' @param op_dir Directory where the files will be saved.
#'
#' @return no return value the objects are saved to the output directory
#' @export
#'
#' @examples
#' \dontrun{
#' save_objects(dds, op_dir = ".")
#' }
save_objects <- function(obj, op_dir = NULL) {
  if(is.null(op_dir)) {
    op_dir = getwd()
  }
  curr_date <- Sys.time() |>
    format("%Y-%M-%d_%H-%m")
  message("Saving the computed objects to: ", op_dir)
  saveRDS(obj, file = file.path(op_dir, paste(class(obj)[1], "-", curr_date, ".RDS", sep = "")))
}
