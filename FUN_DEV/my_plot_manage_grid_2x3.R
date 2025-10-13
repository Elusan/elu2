#' ELU2 management dashboard (2×3) for fitted or managed SPiCT objects
#'
#' @title 2×3 management grid: Biomass, F, Catch, B/Bmsy, F/Fmsy, Kobe
#'
#' @description
#' Builds a compact 2×3 grid of management diagnostics for a **SPiCT** object,
#' accepting **either** a fitted report (`spictcls`) like `m1` **or** the same
#' object after \code{manage()} (i.e., with \code{$man}) like `rman`.
#'
#' **Layout**
#' \itemize{
#'   \item \strong{Row 1}: \code{my_plot_manage_biomass_panel()} \emph{(legend ON)},
#'         \code{my_plot_manage_f_panel()}, \code{my_plot_manage_catch_panel()}
#'   \item \strong{Row 2}: \code{my_plot_manage_bbmsy_panel()},
#'         \code{my_plot_manage_ffmsy_panel()},
#'         \code{my_plot_kobe_all_management_scenario()} \emph{(legends OFF)}
#' }
#'
#' The wrapper preserves each panel’s native behaviour; it only harmonizes
#' legends so that the \emph{first} panel (Biomass) carries the scenario legend.
#'
#' @param obj A fitted \code{spictcls} object (e.g., \code{m1}) or the same
#'   object after \code{manage()} (e.g., \code{rman} with \code{$man} present).
#' @param save Logical; if \code{TRUE}, save a PNG to \code{out_dir}. Default \code{FALSE}.
#' @param out_dir Directory where the PNG is written when \code{save = TRUE}.
#'   Default \code{"FIG/Manage"}.
#' @param file_stub Optional file name (without extension). If \code{NULL}, a
#'   sensible name is generated from whether \code{obj} has \code{$man}.
#' @param width,height,dpi Graphics device size (inches) and resolution for saving.
#'   Defaults: \code{width = 18}, \code{height = 10}, \code{dpi = 300}.
#'
#' @return A patchwork object (the 2×3 grid). If \code{save = TRUE}, the PNG is
#'   written to disk and the patchwork object is still returned (invisibly).
#'
#' @details
#' - This function is a thin orchestrator: each panel is produced by your existing
#'   species functions. It does not alter scales, colors, or data logic.
#' - Only the Biomass panel shows the legend; all other panels hide their legends
#'   to avoid duplication. The Kobe panel is called with legends disabled.
#'
#' @examples
#' \dontrun{
#'   # m1   <- all_models$S1$`S1P.SDM`   # fitted spictcls
#'   # rman <- manage(m1, scenarios = 1:8)
#'
#'   # Works with managed or fitted objects:
#'   # g1 <- my_plot_manage_grid_2x3(m1)
#'   # g2 <- my_plot_manage_grid_2x3(rman, save = TRUE, out_dir = "FIG/Manage")
#'   # print(g1); print(g2)
#' }
#'
#' @seealso
#'   \code{\link{my_plot_manage_biomass_panel}},
#'   \code{\link{my_plot_manage_f_panel}},
#'   \code{\link{my_plot_manage_catch_panel}},
#'   \code{\link{my_plot_manage_bbmsy_panel}},
#'   \code{\link{my_plot_manage_ffmsy_panel}},
#'   \code{\link{my_plot_kobe_all_management_scenario}}
#'
#' @import ggplot2
#' @import patchwork
#' @export
my_plot_manage_grid_2x3 <- function(obj,
                                    save = FALSE,
                                    out_dir = file.path("FIG", "Manage"),
                                    file_stub = NULL,
                                    width = 18, height = 10, dpi = 300) {
  if (!inherits(obj, "spictcls")) {
    stop("`obj` must be a SPiCT report object (class 'spictcls').")
  }
  # Make sure patchwork is available (for safety; we're importing it anyway)
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }

  # ---- Build panels (legend only in the first panel as requested) ----
  p_biomass <- my_plot_manage_biomass_panel(obj, show_legend = TRUE)
  p_f       <- my_plot_manage_f_panel(obj, show_legend = FALSE)
  p_catch   <- my_plot_manage_catch_panel(obj, show_legend = FALSE)

  p_bbmsy   <- my_plot_manage_bbmsy_panel(obj, show_legend = FALSE)
  p_ffmsy   <- my_plot_manage_ffmsy_panel(obj, show_legend = FALSE)

  # Kobe: hide legends here to avoid duplication (legend is on Biomass panel)
  p_kobe <- my_plot_kobe_all_management_scenario(
    obj,
    plot.legend = FALSE,
    man.legend  = FALSE
  )

  # ---- Compose the 2x3 grid ----
  top_row    <- p_biomass | p_f | p_catch
  bottom_row <- p_bbmsy   | p_ffmsy | p_kobe
  grid_2x3   <- (top_row / bottom_row) + patchwork::plot_layout(heights = c(1, 1))

  # ---- Optional save ----
  if (isTRUE(save)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (is.null(file_stub) || !nzchar(file_stub)) {
      type_tag <- if ("man" %in% names(obj)) "managed" else "fitted"
      file_stub <- paste0("manage_grid_2x3_", type_tag)
    }
    out_file <- file.path(out_dir, paste0(file_stub, ".png"))
    ggplot2::ggsave(
      filename = out_file,
      plot     = grid_2x3,
      width    = width,
      height   = height,
      dpi      = dpi
    )
    invisible(grid_2x3)
  } else {
    grid_2x3
  }
}
