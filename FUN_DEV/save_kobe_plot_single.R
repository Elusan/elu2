#' Save a single Kobe plot to disk (robust)
#'
#' @description
#' Builds a Kobe plot using [elu2::kobe_safe()] (which wraps
#' [elu2::kobe_all_in_one_gg()]) and saves it as a PNG.
#' The plot is annotated with a bold top-left ID via [elu2::add_id_label()].
#'
#' @param rep_obj A fitted SPiCT-like object accepted by `kobe_all_in_one_gg()`.
#' @param id_text Character label (e.g., "S1P.SDM") drawn at top-left inside the panel.
#' @param file_path Full output PNG path (directories should exist; use `dir.create(recursive=TRUE)`).
#' @param width,height Numeric, **inches** for PNG device. Default `6 x 5`.
#' @param dpi Numeric PNG resolution. Default `300`.
#' @param ... Passed to [elu2::kobe_safe()] (e.g., `rel.axes`, `logax`, `CI`, etc.).
#'
#' @return Invisibly returns `TRUE` on success, `FALSE` (with a warning) on failure.
#'
#' @import ggplot2
#' @importFrom grDevices png dev.off
#' @export
save_kobe_plot_single <- function(rep_obj,
                                  id_text,
                                  file_path,
                                  width = 6,
                                  height = 5,
                                  dpi = 300,
                                  ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")

  # Build plot using your safe wrapper + add a bold ID in the top-left
  p <- add_id_label(
    kobe_safe(rep_obj, ...),
    id_text = id_text,
    id_size = 4
  )

  # Robust save (1) ggsave then (2) base png fallback
  ok <- try({
    ggplot2::ggsave(
      filename = file_path,
      plot     = p,
      device   = "png",
      width    = width,
      height   = height,
      units    = "in",
      dpi      = dpi,
      bg       = "white"
    )
    TRUE
  }, silent = TRUE)

  if (!inherits(ok, "try-error")) return(invisible(TRUE))

  gr <- try({
    grDevices::png(filename = file_path, width = width, height = height,
                   units = "in", res = dpi, type = "cairo")
    print(p)
    grDevices::dev.off()
    TRUE
  }, silent = TRUE)

  if (inherits(gr, "try-error")) {
    warning("Failed to save plot to: ", file_path)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}
