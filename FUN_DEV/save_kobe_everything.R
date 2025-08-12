#' Save both per-model Kobe plots and the scenario grid
#'
#' @param all_models Named list of scenarios (each a list of fitted models).
#' @param scenario_name Scenario to render (must exist in `names(all_models)`).
#' @param layout "3x2" or "3x1" grid layout for the scenario panel.
#' @param per_model_dir Directory for individual PNGs, e.g., "FIG/per_model".
#' @param grid_dir Directory for per-scenario grid PNGs, e.g., "FIG/KobePhases".
#' @param grid_width,grid_height Grid PNG size in inches (defaults match your build_kobe_grids_scenario()).
#' @param grid_dpi DPI for saved PNGs.
#' @param ... Passed to `kobe_safe()`/`build_kobe_01()`.
#' @return Invisibly, a list with elements `per_model_paths` (named by model)
#'   and `grid_path` for the scenario.
#' @export
save_kobe_everything <- function(all_models,
                                 scenario_name,
                                 layout = "3x2",
                                 per_model_dir = file.path("FIG", "per_model"),
                                 grid_dir      = file.path("FIG", "KobePhases"),
                                 grid_width    = if (layout == "3x2") 15 else 6,
                                 grid_height   = if (layout == "3x2") 6.5 else 10,
                                 grid_dpi      = 300,
                                 ...) {
  if (!scenario_name %in% names(all_models)) {
    stop("Scenario '", scenario_name, "' not found in all_models.")
  }
  sc <- all_models[[scenario_name]]

  # 1) Save per-model PNGs
  scen_model_dir <- file.path(per_model_dir, scenario_name)
  if (!dir.exists(scen_model_dir)) dir.create(scen_model_dir, recursive = TRUE, showWarnings = FALSE)

  pm_paths <- character(0)
  for (nm in names(sc)) {
    p  <- kobe_safe(sc[[nm]], ...)
    fp <- file.path(scen_model_dir, paste0(nm, ".png"))
    .save_plot_safely(p, fp, width = 6, height = 5, dpi = 300)
    pm_paths[nm] <- fp
  }

  # 2) Build and save the scenario grid
  if (!dir.exists(grid_dir)) dir.create(grid_dir, recursive = TRUE, showWarnings = FALSE)
  g  <- build_kobe_01(all_models, scenario_name, layout = layout, ...)
  gp <- file.path(grid_dir, paste0("KobePhases_", scenario_name, "_", layout, ".png"))
  .save_plot_safely(g, gp, width = grid_width, height = grid_height, dpi = grid_dpi)

  invisible(list(per_model_paths = pm_paths, grid_path = gp))
}
