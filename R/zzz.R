#' @useDynLib elu2
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  dll <- system.file("libs", paste0("elu2_tmb", .Platform$dynlib.ext), package = pkgname)
  if (file.exists(dll) && !is.loaded("objective_function")) {
    dyn.load(dll)
  }
}
