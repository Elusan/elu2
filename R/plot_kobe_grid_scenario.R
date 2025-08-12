#' Plot a 2×3 Kobe grid for one scenario (and optionally save)
#'
#' @description
#' Builds a patchwork grid for a single scenario: top row shows SDM models
#' (P, S, F) and bottom row shows GLM models (P, S, F). Each panel is produced
#' via [elu2::kobe_safe()] and tagged with its model ID using
#' [elu2::add_id_label()]. If `save=TRUE`, also saves a PNG to `out_dir`.
#'
#' @param all_models Named list of scenarios. Each scenario is a named list of
#'   fitted model objects with names like `"S1P.SDM"`, `"S1S.GLM"`, etc.
#' @param scenario_name Scenario to render (must exist in `all_models`).
#' @param save Logical; save a PNG of the grid. Default `FALSE`.
#' @param out_dir Directory to save the PNG if `save=TRUE`. Default
#'   `file.path("FIG","KobePhases")`.
#' @param width,height Numeric; inches for the PNG. Default `12 x 5.5`.
#' @param dpi Numeric; dots per inch. Default `300`.
#' @param man.legend Logical; show management-scenario legends in panels.
#'   Default `FALSE` (helps avoid a legend-row changing the layout).
#' @param ... Additional arguments forwarded to [elu2::kobe_safe()]
#'   (e.g., `rel.axes`, `CI`, `logax`, etc.).
#'
#' @return A patchwork object (the composed 2×3 grid). Invisibly returns the
#'   object even if `save=TRUE`.
#'
#' @examples
#' \dontrun{
#' g <- plot_kobe_grid_scenario(all_models, "S1",
#'                              rel.axes = FALSE, CI = 0.95)
#' print(g)
#' }
#' @export
plot_kobe_grid_scenario <- function(all_models, scenario_name,
                                    save = FALSE,
                                    out_dir = file.path("FIG", "KobePhases"),
                                    width = 12, height = 5.5, dpi = 300,
                                    man.legend = FALSE,
                                    ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required.")
  if (!scenario_name %in% names(all_models)) {
    stop("Scenario '", scenario_name, "' not found in all_models.")
  }

  sc <- all_models[[scenario_name]]
  if (!is.list(sc) || !length(sc)) stop("Empty scenario list: ", scenario_name)

  # helpers to pick and tag
  .pick_model <- function(lst, model_letter, data_type) {
    nms <- names(lst)
    idx <- grep(paste0(model_letter, "\\.", data_type, "$"), nms)
    if (length(idx) == 0) return(NULL)
    lst[[idx[1]]]
  }
  get_id <- function(letter, dtype) {
    nm  <- names(sc)
    hit <- nm[grep(paste0(letter, "\\.", dtype, "$"), nm)]
    if (length(hit)) hit[1] else paste0(scenario_name, letter, ".", dtype)
  }

  reps <- list(
    P_SDM = .pick_model(sc, "P", "SDM"),
    S_SDM = .pick_model(sc, "S", "SDM"),
    F_SDM = .pick_model(sc, "F", "SDM"),
    P_GLM = .pick_model(sc, "P", "GLM"),
    S_GLM = .pick_model(sc, "S", "GLM"),
    F_GLM = .pick_model(sc, "F", "GLM")
  )

  # build 6 panels
  p_P_SDM <- elu2::add_id_label(elu2::kobe_safe(reps$P_SDM, man.legend = man.legend, ...),
                                get_id("P","SDM"))
  p_S_SDM <- elu2::add_id_label(elu2::kobe_safe(reps$S_SDM, man.legend = man.legend, ...),
                                get_id("S","SDM"))
  p_F_SDM <- elu2::add_id_label(elu2::kobe_safe(reps$F_SDM, man.legend = man.legend, ...),
                                get_id("F","SDM"))

  p_P_GLM <- elu2::add_id_label(elu2::kobe_safe(reps$P_GLM, man.legend = man.legend, ...),
                                get_id("P","GLM"))
  p_S_GLM <- elu2::add_id_label(elu2::kobe_safe(reps$S_GLM, man.legend = man.legend, ...),
                                get_id("S","GLM"))
  p_F_GLM <- elu2::add_id_label(elu2::kobe_safe(reps$F_GLM, man.legend = man.legend, ...),
                                get_id("F","GLM"))

  # Force a 2 × 3 layout explicitly: row1 = SDM (P,S,F), row2 = GLM (P,S,F)
  design <- "
ABC
DEF
"
  g <- (p_P_SDM + p_S_SDM + p_F_SDM + p_P_GLM + p_S_GLM + p_F_GLM) +
    patchwork::plot_layout(design = design, guides = "keep")

  if (isTRUE(save)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fpath <- file.path(out_dir, paste0("KobeGrid_", scenario_name, "_2x3.png"))
    ok <- try({
      ggplot2::ggsave(filename = fpath, plot = g, device = "png",
                      width = width, height = height, units = "in",
                      dpi = dpi, bg = "white")
      TRUE
    }, silent = TRUE)
    if (inherits(ok, "try-error")) {
      gr <- try({
        grDevices::png(filename = fpath, width = width, height = height,
                       units = "in", res = dpi, type = "cairo")
        print(g)
        grDevices::dev.off()
        TRUE
      }, silent = TRUE)
      if (inherits(gr, "try-error")) warning("Failed to save grid: ", fpath)
    }
  }

  return(invisible(g))
}
