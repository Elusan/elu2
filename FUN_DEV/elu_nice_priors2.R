#' Save Prior vs Posterior Plots for Multiple Models
#'
#' Generates and saves prior vs posterior distribution plots (as PNG files) for each fitted model object in a list.
#'
#' @param models A named list of fitted model objects (e.g., results from \code{fit.elu2()} or similar).
#' @param priors_fun Function to plot priors and posteriors for a single model. Default is \code{priors.elu6}.
#' @param outdir Output directory to save PNG figures. Created if it does not exist. Default: "FIG".
#' @param width Width of the output PNG figure in inches. Default: 8.
#' @param height Height of the output PNG figure in inches. Default: 4.
#' @param dpi Resolution of the PNG figure in dots per inch. Default: 300.
#'
#' @details
#' For each model in \code{models}, this function calls the specified \code{priors_fun} (typically a prior vs posterior plotting function)
#' and saves the resulting plot as a PNG in the output directory.
#'
#' The PNG files are named \code{"priors_<model_name>.png"} according to the names in the input list.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of saving plots.
#'
#' @examples
#' \dontrun{
#' models <- list(
#'   S1F = fit_fox,
#'   S1P = fit_pella,
#'   S1S = fit_schaefer
#' )
#' elu_nice_priors(models)
#' }
#' @export
elu_nice_priors2 <- function(models, priors_fun = priors.elu66, outdir = "FIG", width = 8, height = 4, dpi = 300) {
  if (!dir.exists(outdir)) dir.create(outdir)
  for (nm in names(models)) {
    cat("Saving priors for model:", nm, "\n")
    plot_obj <- priors_fun(models[[nm]], model_id = nm, do.plot = NULL, stamp = NULL)
    ggsave(filename = file.path(outdir, paste0("priors_", nm, ".png")),
           plot = plot_obj, width = width, height = height, dpi = dpi)
  }
}
