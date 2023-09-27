source("renv/activate.R")

### SLUSHY RPROFILE - START ###
## Set slushy global options and environment variables
Sys.setenv(RENV_CONFIG_USE_CACHE = TRUE)
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = TRUE)
## Code to activate renv
local({
  options(renv.config.synchronized.check = FALSE)
    activated_ok <- tryCatch({
       ok <- utils::capture.output({source("renv/activate.R")}, type = "output")
       TRUE
    }, warning = function(w){
       FALSE
  })
  if (!activated_ok){
    lib <- renv::paths$library()
    stop(paste0("Problem activating {renv}. May be due to incorrect version",
                 " of {renv} already installed to project library.",
                 "\nPlease try the following command followed by a session restart:",
                 "\n`renv::remove(\"renv\", library = \"", lib, "\")`"), call. = FALSE)
  }
})
## Code to install correct version of slushy if needed
renv::restore(packages = "slushy", prompt = FALSE)
## Code to check slushy status
local({
  slushy::slushy_status(pkg_deps_ok = TRUE)
})
### SLUSHY RPROFILE - END ###
