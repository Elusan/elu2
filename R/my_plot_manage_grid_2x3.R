#' ELU2 management dashboard (2×3) for fitted or managed SPiCT objects
#'
#' @title 2×3 management grid: Biomass, F, Catch, B/Bmsy, F/Fmsy, Kobe
#'
#' @description
#' Builds a compact 2×3 grid of management diagnostics for a **SPiCT** object,
#' accepting either a fitted report (`spictcls`, e.g. `m1`) or the same object
#' after `manage()` (i.e., with `$man`, e.g. `rman`).
#'
#' **Layout**
#' - **Row 1:** `my_plot_manage_biomass_panel()` (legend ON only if scenarios exist),
#'   `my_plot_manage_f_panel()`, `my_plot_manage_catch_panel()`
#' - **Row 2:** `my_plot_manage_bbmsy_panel()`, `my_plot_manage_ffmsy_panel()`,
#'   `my_plot_kobe_all_management_scenario()` (legend OFF)
#'
#' The wrapper preserves each panel’s native behaviour. The Biomass panel carries
#' the scenario legend when scenarios exist; all other panels hide their legends.
#'
#' @param obj A fitted `spictcls` object (e.g., `m1`) or the same object after
#'   `manage()` (e.g., `rman` with `$man` present).
#' @param save Logical; if `TRUE`, save a PNG to `out_dir`. Default `FALSE`.
#' @param out_dir Directory where the PNG is written when `save = TRUE`.
#'   Default `file.path("FIG","Manage")`.
#' @param file_stub Optional file name (without extension). If `NULL`, a sensible
#'   name is generated from whether `obj` has `$man`.
#' @param width,height,dpi Graphics device size (inches) and resolution for saving.
#'   Defaults: `width = 18`, `height = 10`, `dpi = 300`.
#'
#' @return A patchwork object (the 2×3 grid). If `save = TRUE`, the PNG is
#'   written and the patchwork object is returned invisibly.
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
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }

  # Scenario presence (drives legend visibility on the Biomass panel)
  has_man <- ("man" %in% names(obj)) && length(obj$man) > 0

  # --- Build the six panels ---
  p_biomass <- my_plot_manage_biomass_panel(obj, show_legend = has_man)
  p_f       <- my_plot_manage_f_panel(obj,       show_legend = FALSE)
  p_catch   <- my_plot_manage_catch_panel(obj,   show_legend = FALSE)

  p_bbmsy   <- my_plot_manage_bbmsy_panel(obj,   show_legend = FALSE)
  p_ffmsy   <- my_plot_manage_ffmsy_panel(obj,   show_legend = FALSE)

  # Kobe panel: your updated version supports both roles; keep legends off here
  p_kobe <- my_plot_kobe_all_management_scenario(
    obj,
    plot.legend = FALSE,
    man.legend  = FALSE
  )

  # --- Compose the 2×3 grid ---
  top_row    <- p_biomass | p_f | p_catch
  bottom_row <- p_bbmsy   | p_ffmsy | p_kobe
  grid_2x3   <- (top_row / bottom_row) + patchwork::plot_layout(heights = c(1, 1))

  # --- Optional save ---
  if (isTRUE(save)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (is.null(file_stub) || !nzchar(file_stub)) {
      type_tag <- if (has_man) "managed" else "fitted"
      file_stub <- paste0("manage_grid_2x3_", type_tag)
    }
    out_file <- file.path(out_dir, paste0(file_stub, ".png"))
    ggplot2::ggsave(filename = out_file, plot = grid_2x3,
                    width = width, height = height, dpi = dpi)
    return(invisible(grid_2x3))
  }

  grid_2x3
}
