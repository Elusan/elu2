#' @useDynLib elu2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  dll <- system.file("libs", paste0("elu2", .Platform$dynlib.ext), package = pkgname)
  if (file.exists(dll) && !is.loaded("objective_function")) {
    dyn.load(dll)
  }
}

