elu_nice_priors_999 <- function(models, priors_fun = priors.elu6, outdir = "FIG", width = 8, height = 4, dpi = 300) {
  if (!dir.exists(outdir)) dir.create(outdir)
  for (nm in names(models)) {
    cat("Saving priors for model:", nm, "\n")
    plot_obj <- priors_fun(models[[nm]], model_id = nm, do.plot = NULL, stamp = NULL)
    ggsave(filename = file.path(outdir, paste0("priors_", nm, ".png")),
           plot = plot_obj, width = width, height = height, dpi = dpi)
  }
}
