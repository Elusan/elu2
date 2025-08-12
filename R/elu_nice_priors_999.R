#' Save Prior vs Posterior Plots for Multiple SPiCT Models
#'
#' Generates and saves prior vs posterior distribution plots (as PNG files) for each fitted model object in a list.
#'
#' @param models A named list of fitted SPiCT model objects (e.g., results from \code{fit.elu2()} or similar).
#' @param priors_fun A function used to generate the prior vs posterior plot for a single model.
#'   Default is \code{\link{priors.elu6}}.
#' @param outdir Output directory to save the PNG figures. If it does not exist, it will be created.
#'   Default: \code{"FIG"}.
#' @param width Width of each PNG figure in inches. Default: \code{8}.
#' @param height Height of each PNG figure in inches. Default: \code{4}.
#' @param dpi Resolution of the output PNG figure in dots per inch. Default: \code{300}.
#'
#' @details
#' For each model in the input list, the function generates a prior–posterior comparison plot using
#' the specified \code{priors_fun} and saves it as a PNG file named \code{priors_<model>.png}
#' in the specified output directory.
#'
#' @return Invisibly returns \code{NULL}. Side effect: saves PNG files to disk.
#'
#' @examples
#' \dontrun{
#' models <- list(S1P = fit1, S1F = fit2, S1S = fit3)
#' elu_nice_priors_999(models)
#' }
#'
#' @export
elu_nice_priors_999 <- function(models, priors_fun = priors.elu6, outdir = "FIG/Elu_Priors", width = 8, height = 4, dpi = 300) {
  if (!dir.exists(outdir)) dir.create(outdir)
  for (nm in names(models)) {
    cat("Saving priors for model:", nm, "\n")
    plot_obj <- priors_fun(models[[nm]], model_id = nm, do.plot = NULL, stamp = NULL)
    ggsave(filename = file.path(outdir, paste0("priors_", nm, ".png")),
           plot = plot_obj, width = width, height = height, dpi = dpi)
  }
}
