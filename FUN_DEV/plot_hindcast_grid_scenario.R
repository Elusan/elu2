#' Hindcast grid (2×3) for one scenario
#'
#' Builds a 2×3 grid: columns = **Pella**, **Schaefer**, **Fox**;
#' rows = **SDM** (top), **GLM** (bottom). Panel labels are like `"S1P.SDM"`.
#' Missing models or plotting errors produce placeholders.
#'
#' @param all_models Named list-of-lists by scenario, e.g. `all_models$S1$S1P.SDM`.
#' @param scenario Name of the scenario to plot, e.g. `"S1"`.
#' @param npeels,peel.dtc,mc.cores Hindcast controls (passed to the internal panel constructor).
#' @param CI,verbose,add.mase Passed through to \code{plotspict.hindcast_elu2_gg_exact_FOR_GRIDS()}.
#'   (Shown for provenance only; this symbol is not linked to avoid roxygen link warnings.)
#' @param show_legends Logical; keep per-panel legends (default `TRUE`). Set `FALSE` to declutter.
#' @param save Logical; save PNG to `out_dir`. Default `TRUE`.
#' @param out_dir Output directory (created if missing). Default `file.path("FIG", "Hindcast")`.
#' @param filename Optional filename; default `"HindcastGrid_<scenario>_2x3.png"`.
#' @param width,height,dpi PNG device settings (inches, inches, dpi).
#'
#' @return Invisibly, a list with elements:
#'   \itemize{
#'     \item \code{patch}: the patchwork object
#'     \item \code{file}: the saved file path (or \code{NA} if \code{save = FALSE})
#'   }
#'
#' @seealso \link[=plot_hindcast_grid_many]{plot_hindcast_grid_many()} for batch rendering
#'   across multiple scenarios.
#'
#' @examples
#' \dontrun{
#' # Suppose all_models is a named list-of-lists as described:
#' res <- plot_hindcast_grid_scenario(all_models, scenario = "S1",
#'                                    npeels = 3, peel.dtc = FALSE,
#'                                    out_dir = file.path("FIG", "Hindcast"),
#'                                    width = 12, height = 5, dpi = 300)
#' res$file  # saved path
#' }
#'
#' @export
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggsave
plot_hindcast_grid_scenario <- function(
    all_models, scenario,
    npeels = 2, peel.dtc = FALSE, mc.cores = 1,
    CI = 0.95, verbose = FALSE, add.mase = TRUE,
    show_legends = TRUE,
    save = TRUE, out_dir = file.path("FIG", "Hindcast"),
    filename = NULL, width = 12, height = 5, dpi = 300
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install it with install.packages('patchwork').")
  }
  if (!scenario %in% names(all_models)) {
    stop("Scenario '", scenario, "' not found in all_models.")
  }
  mods <- all_models[[scenario]]

  # Expected keys by design (columns = P,S,F ; rows = SDM, GLM)
  keys <- c(
    paste0(scenario, c("P", "S", "F"), ".SDM"),
    paste0(scenario, c("P", "S", "F"), ".GLM")
  )

  # Helper: robust panel maker (placeholder on error or missing fit)
  make_safe_panel <- function(k) {
    if (!k %in% names(mods)) {
      message("[", k, "] not found; placing placeholder panel.")
      return(make_hindcast_panel(fit = NULL, panel_label = k))
    }
    fit <- mods[[k]]
    msg <- paste0("[", k, "]")
    p <- try({
      message(msg, " running hindcast & composing panel...")
      make_hindcast_panel(
        fit, panel_label = k,
        npeels = npeels, peel.dtc = peel.dtc, mc.cores = mc.cores,
        CI = CI, verbose = verbose, add.mase = add.mase,
        show_legend = show_legends
      )
    }, silent = TRUE)
    if (inherits(p, "try-error")) {
      message(msg, " plot failed; showing placeholder panel.")
      make_hindcast_panel(fit = NULL, panel_label = k)
    } else {
      p
    }
  }

  panels <- lapply(keys, make_safe_panel)

  patch <- patchwork::wrap_plots(panels, nrow = 2, byrow = TRUE)

  outfile <- NA_character_
  if (isTRUE(save)) {
    if (is.null(filename)) {
      filename <- sprintf("HindcastGrid_%s_2x3.png", scenario)
    }
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    outfile <- file.path(out_dir, filename)
    ggplot2::ggsave(outfile, patch, width = width, height = height, dpi = dpi)
    message("Saved: ", outfile)
  }

  invisible(list(patch = patch, file = outfile))
}


#' Hindcast grid for many scenarios (batch)
#'
#' Iterates scenarios in \code{scenario_names} (default: all found in \code{all_models})
#' and writes one 2×3 grid PNG per scenario using
#' \link[=plot_hindcast_grid_scenario]{plot_hindcast_grid_scenario()}.
#'
#' @param all_models See \link[=plot_hindcast_grid_scenario]{plot_hindcast_grid_scenario()}.
#' @param scenario_names Character vector of scenarios to process; default \code{names(all_models)}.
#' @param ... Passed through to \link[=plot_hindcast_grid_scenario]{plot_hindcast_grid_scenario()}.
#'
#' @return Invisibly, a named list of results (one per scenario), each as returned by
#'   \link[=plot_hindcast_grid_scenario]{plot_hindcast_grid_scenario()}.
#'
#' @examples
#' \dontrun{
#' # Render all scenarios:
#' out <- plot_hindcast_grid_many(all_models,
#'                                out_dir = file.path("FIG", "Hindcast"),
#'                                width = 12, height = 5, dpi = 300,
#'                                show_legends = FALSE)
#' }
#'
#' @export
plot_hindcast_grid_many <- function(all_models, scenario_names = names(all_models), ...) {
  stopifnot(is.list(all_models))
  out <- list()
  for (sc in scenario_names) {
    message("== Scenario ", sc, " ==")
    out[[sc]] <- plot_hindcast_grid_scenario(all_models, scenario = sc, ...)
  }
  invisible(out)
}
