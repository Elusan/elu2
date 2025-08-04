#' Save Prior vs Posterior Plots for Multiple Models, Plus Scenario Grids
#'
#' For each model in \code{models}, saves an individual prior–posterior plot.
#' Then, for each scenario (e.g. "S1"), collects its three submodels
#' ("S1P", "S1F", "S1S"), stacks them as rows, and saves a combined grid.
#'
#' @param models Named list of fitted model objects. Names must end in "P", "F", or "S",
#'   prefixed by scenario tag (e.g. "S1P", "S1F", "S1S").
#' @param priors_fun Function to plot priors vs posterior for one model.
#'   Signature: \code{priors_fun(model_obj, model_id, ...)} returning a ggplot.
#'   Default: \code{priors.elu6}.
#' @param outdir Output directory for PNGs. Default: `"FIG"`.
#' @param width Width (in inches) of each individual plot. Default: 8.
#' @param height Height (in inches) of each individual plot. Default: 4.
#' @param dpi Resolution (dpi) for saved PNGs. Default: 300.
#'
#' @details
#' - Individual plots are named: `priors_<model_name>.png`
#' - Scenario grids are named: `priors_<scenario>_grid.png`,
#'   where `<scenario>` is the prefix (e.g. "S1").
#' - Uses patchwork to stack the three ggplots in row order: Pella, Fox, Schaefer.
#'
#' @return Invisibly \code{NULL}. Called for side effects.
#'
#' @importFrom patchwork wrap_plots
#' @export
elu_nice_priors10 <- function(models,
                            priors_fun = priors.elu6,
                            outdir     = "FIG",
                            width      = 8,
                            height     = 4,
                            dpi        = 300) {
  # ensure output dir
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # 1) Generate and save each individual plot, keep in a list
  plot_list <- list()
  for (nm in names(models)) {
    message("Saving priors for model: ", nm)
    p <- priors_fun(models[[nm]], model_id = nm)
    ggsave(
      filename = file.path(outdir, paste0("priors_", nm, ".png")),
      plot     = p,
      width    = width,
      height   = height,
      dpi      = dpi
    )
    plot_list[[nm]] <- p
  }

  # 2) Identify unique scenarios by stripping trailing model letter
  scenarios <- unique(sub("(P|F|S)$", "", names(models)))

  # 3) For each scenario, stack its P, F, S plots
  for (scn in scenarios) {
    # define the three expected submodel names in order
    subs <- paste0(scn, c("P", "F", "S"))
    subs_present <- subs[subs %in% names(plot_list)]
    if (length(subs_present) == 0) {
      warning("No models found for scenario: ", scn)
      next
    }

    # collect the ggplots in the correct row order
    grid_plots <- plot_list[subs_present]

    # stack vertically: P / F / S  (or however many are present)
    combined <- patchwork::wrap_plots(grid_plots, ncol = 1)

    # save the grid
    grid_height <- height * length(grid_plots)
    out_file <- file.path(outdir, paste0("priors_", scn, "_grid.png"))
    message("Saving scenario grid for: ", scn, " -> ", out_file)
    ggsave(
      filename = out_file,
      plot     = combined,
      width    = width,
      height   = grid_height,
      dpi      = dpi,
      limitsize = FALSE
    )
  }

  invisible(NULL)
}
