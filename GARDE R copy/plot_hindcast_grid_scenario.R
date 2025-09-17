# R/plot_hindcast_grid_scenario.R

#' Hindcast grid (2×3) for one scenario
#'
#' Builds a 2×3 grid: columns = **Pella**, **Schaefer**, **Fox**;
#' rows = **SDM** (top), **GLM** (bottom). Panel labels are like `"S1P.SDM"`.
#' Missing models or plotting errors produce placeholders.
#'
#' @param all_models Named list-of-lists by scenario, e.g. `all_models$S1$S1P.SDM`.
#' @param scenario Name of the scenario to plot, e.g. `"S1"`.
#' @param npeels,peel.dtc,mc.cores Hindcast controls.
#' @param CI,verbose,add.mase Passed to [plotspict.hindcast_elu2_gg_exact_FOR_GRIDS()].
#' @param show_legends Logical; keep per-panel legends (default `TRUE`). Set `FALSE` to declutter.
#' @param save Logical; save PNG to `out_dir`. Default `TRUE`.
#' @param out_dir Output directory (created if missing). Default `file.path("FIG", "Hindcast")`.
#' @param filename Optional filename; default `"HindcastGrid_<scenario>_2x3.png"`.
#' @param width,height,dpi PNG device settings (inches, inches, dpi).
#'
#' @return Invisibly, a list with `patch` (patchwork object) and `file` (path or `NA`).
#'
#' @export
#' @importFrom patchwork wrap_plots
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
  mods <- all_models[`scenario`]

  # expected keys by design (columns = P,S,F ; rows = SDM, GLM)
  keys <- c(
    paste0(scenario, c("P", "S", "F"), ".SDM"),
    paste0(scenario, c("P", "S", "F"), ".GLM")
  )

  panels <- lapply(keys, function(k) {
    if (!k %in% names(mods)) {
      message("[", k, "] not found; placing placeholder panel.")
      return(make_hindcast_panel(fit = NULL, panel_label = k))  # will error -> placeholder
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
  })

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
