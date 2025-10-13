#' Hindcast grid for many scenarios (batch)
#'
#' Iterates scenarios in `scenario_names` (default: all found in `all_models`)
#' and writes one 2×3 grid PNG per scenario using
#' \link[=plot_hindcast_grid_scenario]{plot_hindcast_grid_scenario()}.
#'
#' @param all_models See \link[=plot_hindcast_grid_scenario]{plot_hindcast_grid_scenario()}.
#' @param scenario_names Character vector of scenarios to process; default `names(all_models)`.
#' @param ... Passed through to \link[=plot_hindcast_grid_scenario]{plot_hindcast_grid_scenario()}.
#'
#' @return Invisibly, a named list of results (one per scenario).
#'
#' @seealso
#' \link[=plot_hindcast_grid_scenario]{plot_hindcast_grid_scenario()}.
#'
#' @export
plot_hindcast_grid_many <- function(all_models, scenario_names = names(all_models), ...) {
  out <- list()
  for (sc in scenario_names) {
    message("== Scenario ", sc, " ==")
    out[[sc]] <- plot_hindcast_grid_scenario(all_models, scenario = sc, ...)
  }
  invisible(out)
}
