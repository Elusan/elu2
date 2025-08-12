#' Kobe plotting core
#'
#' Core helpers that sit around \code{kobe_all_in_one_gg()}:
#' - \code{kobe_safe()}: robust single-panel builder that never errors
#' - \code{add_id_label()}: stamp a bold panel ID at top-left
#' - \code{.pick_model()}: internal selector for P/S/F × SDM/GLM
#' - \code{.save_plot_safely()}: internal PNG saver with fallback
#'
#' @name kobe_core
#' @keywords internal
NULL


#' Safely build a single Kobe plot (robust to failures)
#' @rdname kobe_core
#' @import ggplot2
#' @export
kobe_safe <- function(rep_obj,
                      logax = FALSE, plot.legend = TRUE, man.legend = TRUE,
                      ext = TRUE, rel.axes = FALSE,
                      xlim = NULL, ylim = NULL,
                      labpos = c(1, 1), xlabel = NULL, stamp = NULL,
                      verbose = TRUE, CI = 0.95) {
  res <- try(
    {
      # delegate the real work to your ggplot builder
      p <- kobe_all_in_one_gg(
        rep = rep_obj,
        logax = logax, plot.legend = plot.legend, man.legend = man.legend,
        ext = ext, rel.axes = rel.axes,
        xlim = xlim, ylim = ylim,
        labpos = labpos, xlabel = xlabel, stamp = stamp,
        verbose = verbose, CI = CI
      )
      p
    },
    silent = TRUE
  )

  # If anything failed, return a minimal but valid ggplot
  if (inherits(res, "try-error") || is.null(res)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
    p_err <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0, y = 0,
        label = "Plot error",
        size = 4, colour = "red", fontface = "bold"
      ) +
      ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
      ggplot2::theme_void()
    return(p_err)
  }
  res
}


#' Add a model ID label inside a plot (top-left)
#' @rdname kobe_core
#' @import ggplot2
#' @export
add_id_label <- function(p, id_text, id_size = 4) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  p +
    ggplot2::annotate(
      "text",
      x = -Inf, y = Inf, label = id_text,
      hjust = -0.2, vjust = 2,     # small, empirically chosen nudges
      fontface = "bold", size = id_size
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 6))
}


#' @rdname kobe_core
#' @keywords internal
#' @noRd
.pick_model <- function(lst, model_letter, data_type) {
  # Example names like "S1P.GLM", "S1S.SDM", etc.
  # We find the FIRST name in `lst` that ends with "<letter>.<data_type>"
  # and return that element; otherwise NULL.
  if (!length(lst)) return(NULL)
  nms <- names(lst)
  idx <- grep(paste0(model_letter, "\\.", data_type, "$"), nms)
  if (length(idx) == 0) return(NULL)
  lst[[idx[1]]]
}


#' @rdname kobe_core
#' @keywords internal
#' @noRd
#' @importFrom grDevices png dev.off
.save_plot_safely <- function(plot_obj, file_path, width, height, dpi) {
  ok <- try({
    ggplot2::ggsave(
      filename = file_path,
      plot     = plot_obj,
      device   = "png",
      width    = width,
      height   = height,
      units    = "in",
      dpi      = dpi,
      bg       = "white"  # helps with Beamer/Word transparency issues
    )
    TRUE
  }, silent = TRUE)

  if (!inherits(ok, "try-error")) return(invisible(TRUE))

  gr <- try({
    grDevices::png(
      filename = file_path,
      width    = width,
      height   = height,
      units    = "in",
      res      = dpi,
      type     = "cairo"
    )
    print(plot_obj)
    grDevices::dev.off()
    TRUE
  }, silent = TRUE)

  if (inherits(gr, "try-error")) {
    warning("Failed to save plot to: ", file_path)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}
