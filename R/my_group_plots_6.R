#' ELU2: 2×3 summary grid from existing panels (ggplot2)
#'
#' @description
#' Builds a 2×3 grid from your six existing ELU2/SPiCT ggplot panels and prints it.
#' It accepts either a fitted model (e.g., `m1`) or a management object
#' produced by `manage()` (e.g., `rman`) and **does not** change any of the
#' internal settings of your panel functions—those are called as-is.
#'
#' Panels (left→right, top→bottom):
#' 1) `my_plot_manage_biomass_panel()`
#' 2) `my_plot_manage_f_panel()`
#' 3) `my_plot_manage_catch_panel()`
#' 4) `my_plot_manage_bbmsy_panel()`
#' 5) `my_plot_manage_ffmsy_panel()`
#' 6) `my_plot_kobe_all_management_scenario()`
#'
#' @param x A fitted SPiCT report object (like `m1`) or a management object
#'   returned by `manage()` (like `rman`). The same object will be passed to
#'   all six panel functions.
#' @param save Logical; if `TRUE`, saves the composed grid to `out_dir/filename`.
#'   Default `FALSE`.
#' @param out_dir Directory where the image is written when `save = TRUE`.
#'   Default `"FIG"`.
#' @param filename Output filename when `save = TRUE`. If `NULL`, a name is
#'   generated automatically using the object class and current time.
#' @param width,height,dpi Device width/height (inches) and resolution for saving.
#'   Defaults: `width = 20`, `height = 9`, `dpi = 300`.
#' @param show_legend_biomass Show legend in the biomass panel?
#'   Default `TRUE` (matches your example call).
#' @param show_legend_f Show legend in the F panel? Default `FALSE`.
#' @param show_legend_catch Show legend in the catch panel? Default `FALSE`.
#' @param show_legend_bbmsy Show legend in the B/Bmsy panel? Default `TRUE`.
#' @param show_legend_ffmsy Show legend in the F/Fmsy panel? Default `FALSE`.
#' @param guides_collect If `TRUE`, use patchwork legend collection.
#'   Default `FALSE` to avoid altering your panel legend behavior.
#'
#' @return Invisibly returns the patchwork object (and, if `save = TRUE`, the file path as an attribute `"file"`).
#'
#' @examples
#' \dontrun{
#' # Model object:
#' g1 <- my_group_plots_6(m1, save = TRUE, out_dir = "FIG/KobePhasesNEW")
#'
#' # Management object:
#' g2 <- my_group_plots_6(rman)
#' }
#'
#' @export
my_group_plots_6 <- function(
    x,
    save = FALSE,
    out_dir = "FIG",
    filename = NULL,
    width = 20, height = 9, dpi = 300,
    show_legend_biomass = TRUE,
    show_legend_f       = FALSE,
    show_legend_catch   = FALSE,
    show_legend_bbmsy   = TRUE,
    show_legend_ffmsy   = FALSE,
    guides_collect      = FALSE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }

  # --- Build each panel by calling your existing functions WITHOUT changing them ---
  # Use your example legend choices as the defaults but allow override above.
  p_biomass <- my_plot_manage_biomass_panel(x, show_legend = show_legend_biomass)
  p_f       <- my_plot_manage_f_panel      (x, show_legend = show_legend_f)
  p_catch   <- my_plot_manage_catch_panel  (x, show_legend = show_legend_catch)
  p_bbmsy   <- my_plot_manage_bbmsy_panel  (x, show_legend = show_legend_bbmsy)
  p_ffmsy   <- my_plot_manage_ffmsy_panel  (x, show_legend = show_legend_ffmsy)
  p_kobe    <- my_plot_kobe_all_management_scenario(x)

  # --- Compose 2x3 grid (left→right, top→bottom) ---
  # We do not force any theme/legend modifications to “keep same settings”.
  g <- (p_biomass | p_f | p_catch) /
    (p_bbmsy   | p_ffmsy | p_kobe)

  # Optionally collect legends (kept OFF by default to avoid changing layout)
  if (guides_collect) {
    g <- g + patchwork::plot_layout(guides = "collect")
  }

  # Print for side-effect (like plot.spictcls)
  print(g)

  # Optional save
  if (isTRUE(save)) {
    if (!dir.exists(out_dir)) {
      ok <- try(dir.create(out_dir, recursive = TRUE), silent = TRUE)
      if (inherits(ok, "try-error") || !dir.exists(out_dir)) {
        warning("Could not create out_dir = '", out_dir, "'. Skipping save.")
      }
    }
    if (is.null(filename)) {
      # Build a friendly filename from class/time
      obj_class <- try(class(x)[1], silent = TRUE)
      if (inherits(obj_class, "try-error") || !nzchar(obj_class)) obj_class <- "elu2"
      ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
      filename <- paste0("my_group_plots_6_", obj_class, "_", ts, ".png")
    }
    outfile <- file.path(out_dir, filename)
    ggplot2::ggsave(
      filename = outfile,
      plot     = g,
      width    = width,
      height   = height,
      dpi      = dpi,
      limitsize = FALSE
    )
    attr(g, "file") <- outfile
    message("Saved: ", outfile)
  }

  invisible(g)
}
