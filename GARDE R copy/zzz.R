#' @useDynLib elu2, .registration = TRUE
NULL

.onLoad <- function(libname, pkgname) {
  # Bind exported datasets into the package namespace so `elu2::name` works.
  ds <- c(
    "Data_glm_full","Data_glm_short","Data_sdm",
    "Delta_GLM_list","Long_delta_GLM_list","Long_sdmTMB_list",
    "My3_indices_list","inp1.glm","inp1.sdmTMB",
    "long_series_list","scenario_definitions","sdmTMB_inp","short_sdmTMB_list"
  )
  ns <- asNamespace(pkgname)
  for (nm in ds) {
    try(utils::data(list = nm, package = pkgname, envir = ns), silent = TRUE)
  }
}

.onUnload <- function(libpath) {
  # Cleanly unload DLL (since we used @useDynLib with registration).
  library.dynam.unload("elu2", libpath)
}
