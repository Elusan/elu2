addto <- function(fname) {
  from <- file.path("FUN_DEV", fname)
  to   <- file.path("R", fname)
  if (!file.exists(from)) stop("Function file does not exist in FUN_DEV/")
  file.copy(from, to, overwrite = TRUE)
  message(sprintf("Function '%s' copied to R/ folder.", fname))
  message("Remember to add roxygen2 documentation if not already present!")
  message("Then run devtools::document() and devtools::install()")
}
