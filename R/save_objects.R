save_objects <- function(obj) {
  curr_date <- Sys.time() |>
    format("%Y-%M-%d_%H-%m")
  saveRDS(obj, file = paste(class(obj)[1], "-", curr_date, ".RDS", sep = ""))
}
