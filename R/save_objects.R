save_objects <- function(obj) {
  curr_date <- Sys.time() |>
    format("%Y-%M-%d_%H-%m")
  qs::qsave(obj, file = paste(class(obj)[1], "-", curr_date, ".qs", sep = ""))
}
