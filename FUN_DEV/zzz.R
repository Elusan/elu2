#' @useDynLib ELU, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


.onLoad <- function(libname, pkgname) {
  library(TMB)
  dll <- system.file("libs", paste0("elu2_tmb", .Platform$dynlib.ext), package = pkgname)
  if (file.exists(dll)) dyn.load(dll)
}
