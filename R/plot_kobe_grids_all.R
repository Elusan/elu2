#' Plot/save Kobe grids for many scenarios in one go
#'
#' @description
#' Iterates over scenarios and calls [elu2::plot_kobe_grid_scenario()] for each.
#' By default it processes **all** scenarios found in `all_models`, and will save
#' one PNG per scenario to `out_dir`. Returns the patchwork objects (one per
#' scenario) invisibly; if `save = TRUE`, also returns the file paths.
#'
#' @param all_models Named list of scenarios. Each scenario is a named list of
#'   fitted model objects accepted by [elu2::kobe_all_in_one_gg()].
#' @param scenario_names Character vector of scenario names to process.
#'   Default: `names(all_models)` (all scenarios).
#' @param save Logical; save each grid PNG to `out_dir`. Default `TRUE`.
#' @param out_dir Directory where PNGs are written when `save = TRUE`.
#'   Default `file.path("FIG", "KobePhases")`.
#' @param width,height Numeric; inches for each grid PNG. Default `12 x 5.5`.
#' @param dpi Numeric; dots per inch for PNGs. Default `300`.
#' @param man.legend Logical; show management-scenario legends inside panels.
#'   Default `FALSE` to avoid an extra legend row changing the layout.
#' @param verbose Logical; emit progress messages. Default `TRUE`.
#' @param ... Additional arguments forwarded to [elu2::plot_kobe_grid_scenario()]
#'   and ultimately to [elu2::kobe_safe()] (e.g., `rel.axes`, `CI`, `logax`).
#'
#' @return Invisibly returns a list with:
#'   * `grids` — named list of patchwork objects (one per scenario)
#'   * `paths` — named character vector of file paths (only when `save = TRUE`)
#'
#' @examples
#' \dontrun{
#' # Build & save grids for all scenarios
#' plot_kobe_grids_all(all_models,
#'                     save = TRUE, out_dir = "FIG/KobePhases",
#'                     width = 12, height = 5.5, dpi = 300,
#'                     man.legend = FALSE,
#'                     rel.axes = FALSE, CI = 0.95)
#'
#' # Only build (no files written), and for a subset:
#' res <- plot_kobe_grids_all(all_models, scenario_names = c("S1","S3"),
#'                            save = FALSE, rel.axes = FALSE, CI = 0.95)
#' print(res$grids$S1)
#' }
#' @export
plot_kobe_grids_all <- function(all_models,
                                scenario_names = names(all_models),
                                save = TRUE,
                                out_dir = file.path("FIG", "KobePhases"),
                                width = 12, height = 5.5, dpi = 300,
                                man.legend = FALSE,
                                verbose = TRUE,
                                ...) {
  if (is.null(scenario_names) || !length(scenario_names)) {
    stop("No scenarios to process: 'scenario_names' is empty.")
  }
  missing_sc <- setdiff(scenario_names, names(all_models))
  if (length(missing_sc)) {
    stop("These scenarios are not in 'all_models': ",
         paste(missing_sc, collapse = ", "))
  }

  if (isTRUE(save) && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  grids <- vector("list", length(scenario_names))
  names(grids) <- scenario_names
  paths <- if (isTRUE(save)) character(length(scenario_names)) else NULL
  if (!is.null(paths)) names(paths) <- scenario_names

  i <- 1L
  while (i <= length(scenario_names)) {
    scen <- scenario_names[i]
    if (isTRUE(verbose)) message("Building grid for scenario: ", scen)

    g <- elu2::plot_kobe_grid_scenario(all_models, scen,
                                       save     = save,
                                       out_dir  = out_dir,
                                       width    = width,
                                       height   = height,
                                       dpi      = dpi,
                                       man.legend = man.legend,
                                       ...)

    grids[[i]] <- g
    if (isTRUE(save)) {
      # Keep filename aligned with plot_kobe_grid_scenario()
      paths[i] <- file.path(out_dir, paste0("KobeGrid_", scen, "_2x3.png"))
    }
    i <- i + 1L
  }

  invisible(list(grids = grids, paths = paths))
}
