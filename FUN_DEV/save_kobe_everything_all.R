#' Save Kobe outputs for many scenarios in one go
#'
#' @description
#' Convenience wrapper that iterates over scenarios and calls
#' [elu2::save_kobe_everything()] for each. By default it processes
#' **all** scenarios in `all_models`. It saves per-model PNGs under
#' `FIG/per_model/<scenario>/` and one grid PNG per scenario under
#' `FIG/KobePhases/`.
#'
#' @param all_models Named list of scenarios. Each scenario is a list of fitted
#'   models accepted by [elu2::kobe_all_in_one_gg()].
#' @param scenario_names Character vector of scenario names to process.
#'   Default: `names(all_models)` (all scenarios).
#' @param layout `"3x2"` (default) or `"3x1"` grid layout for the scenario panel.
#' @param per_model_dir Directory for individual PNGs (default `"FIG/per_model"`).
#' @param grid_dir Directory for per-scenario grid PNGs (default `"FIG/KobePhases"`).
#' @param grid_width,grid_height Grid PNG size in inches (defaults match
#'   `build_kobe_grids_scenario()` for the chosen `layout`).
#' @param grid_dpi DPI for grid PNGs (default `300`).
#' @param ... Additional arguments forwarded to [elu2::kobe_safe()] /
#'   [elu2::build_kobe_01()] via [elu2::save_kobe_everything()] (e.g., `rel.axes`, `CI`, etc.).
#'
#' @return Invisibly, a named list where each element (per scenario) contains:
#'   `per_model_paths` (named character vector of PNG paths) and `grid_path` (PNG path).
#'
#' @examples
#' \dontrun{
#' # Process *all* scenarios
#' save_kobe_everything_all(all_models, layout = "3x2", rel.axes = FALSE, CI = 0.95)
#'
#' # Or a subset
#' save_kobe_everything_all(all_models, scenario_names = c("S1","S3","S6"),
#'                          layout = "3x1", rel.axes = FALSE, CI = 0.95)
#' }
#' @export
save_kobe_everything_all <- function(all_models,
                                     scenario_names = names(all_models),
                                     layout = "3x2",
                                     per_model_dir = file.path("FIG", "per_model"),
                                     grid_dir      = file.path("FIG", "KobePhases"),
                                     grid_width    = if (layout == "3x2") 15 else 6,
                                     grid_height   = if (layout == "3x2") 6.5 else 10,
                                     grid_dpi      = 300,
                                     ...) {
  if (is.null(scenario_names) || !length(scenario_names)) {
    stop("No scenarios to process: 'scenario_names' is empty.")
  }
  missing_sc <- setdiff(scenario_names, names(all_models))
  if (length(missing_sc)) {
    stop("These scenarios are not in 'all_models': ",
         paste(missing_sc, collapse = ", "))
  }

  # Ensure top-level dirs exist up-front (per-scenario subdirs created downstream)
  if (!dir.exists(per_model_dir)) dir.create(per_model_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(grid_dir))      dir.create(grid_dir,      recursive = TRUE, showWarnings = FALSE)

  out <- vector("list", length(scenario_names))
  names(out) <- scenario_names

  i <- 1L
  while (i <= length(scenario_names)) {
    scen <- scenario_names[i]
    # Delegate to your existing single-scenario function
    res <- save_kobe_everything(all_models,
                                scenario_name = scen,
                                layout        = layout,
                                per_model_dir = per_model_dir,
                                grid_dir      = grid_dir,
                                grid_width    = grid_width,
                                grid_height   = grid_height,
                                grid_dpi      = grid_dpi,
                                ...)
    out[[i]] <- res
    i <- i + 1L
  }

  invisible(out)
}
