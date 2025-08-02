#' Save Prior vs Posterior Plots for All Scenarios and Models
#'
#' Handles nested structure from `models_all` and saves PNGs of prior-posterior plots for each fitted model.
#'
#' @param models_all A nested list: models_all[[scenario]][[model]] = list(result = ..., error = ...)
#' @param priors_fun Function to generate prior vs posterior plots. Default: priors.elu6
#' @param outdir Directory to save figures. Default: "FIG"
#' @param width Width in inches. Default: 8
#' @param height Height in inches. Default: 4
#' @param dpi Image resolution. Default: 300
#'
#' @return Invisibly returns NULL
#' @export
elu_nice_priors_all <- function(models_all, priors_fun = priors.elu6,
                                outdir = "FIG/Priors", width = 8, height = 6, dpi = 300) {
  if (!dir.exists(outdir)) dir.create(outdir)

  for (scenario in names(models_all)) {
    for (model_tag in names(models_all[[scenario]])) {
      model_obj <- models_all[[scenario]][[model_tag]]
      if (!is.null(model_obj$error)) {
        message(sprintf("Skipping %s (error in fit)", model_tag))
        next
      }
      if (is.null(model_obj$result)) {
        message(sprintf("Skipping %s (no result)", model_tag))
        next
      }

      cat("Saving priors for model:", model_tag, "\n")
      plot_obj <- priors_fun(model_obj$result, model_id = model_tag)
      ggsave(filename = file.path(outdir, paste0("priors_", model_tag, ".png")),
             plot = plot_obj, width = width, height = height, dpi = dpi)
    }
  }

  invisible(NULL)
}
