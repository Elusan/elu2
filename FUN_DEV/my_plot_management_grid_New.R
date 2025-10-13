#' 2×3 management grid (biomass, F, catch, B/Bmsy, F/Fmsy, Kobe)
#'
#' @description
#' Builds a 2×3 grid of management panels for a fitted **SPiCT** report object:
#' first row shows absolute **Biomass**, **Fishing mortality F[t]**, and **Catch**;
#' second row shows **B/Bmsy**, **F/Fmsy**, and the **Kobe phase plot**.
#'
#' The function delegates to your six panel builders and preserves each panel’s
#' own defaults and appearance. Legends are not merged or re-styled unless
#' explicitly requested. When \code{engine = "cowplot"} and management scenarios
#' are present (\code{rep$man}), the scenario legend is shown **only** on the
#' Absolute Biomass panel by default and suppressed on the other five panels;
#' the Kobe plot keeps its own legend behavior.
#'
#' @param rep A fitted SPiCT report object (\code{spictcls}). May include
#'   management scenarios in \code{rep$man}.
#' @param title Optional title drawn above the whole grid (engine dependent).
#' @param widths Numeric vector of length 3 with relative column widths. Default \code{c(1,1,1)}.
#' @param heights Numeric vector of length 2 with relative row heights. Default \code{c(1,1)}.
#' @param engine Which layout backend to use: \code{"cowplot"}, \code{"patchwork"}, or \code{"gridExtra"}.
#' @param patchwork_guides For \code{engine = "patchwork"}, whether to \code{"keep"}
#'   per-panel guides or \code{"collect"} them. Default \code{"keep"}.
#' @param ... Additional arguments forwarded to the underlying panel functions.
#'   Only arguments matching a panel’s formal parameters are passed through to that panel.
#'
#' @details
#' Panels used (looked up at runtime; errors render a labeled placeholder panel):
#' \itemize{
#'   \item \code{my_plot_manage_biomass_panel(rep_man = rep, ...)}
#'   \item \code{my_plot_manage_f_panel(rep_man = rep, ...)}
#'   \item \code{my_plot_manage_catch_panel(rep_man = rep, ...)}
#'   \item \code{my_plot_manage_bbmsy_panel(rep_man = rep, ...)}
#'   \item \code{my_plot_manage_ffmsy_panel(rep_man = rep, ...)}
#'   \item \code{my_plot_kobe_all_management_scenario(rep = rep, ...)}
#' }
#'
#' @return
#' A plot object suitable for printing:
#' \itemize{
#'   \item For \code{engine = "cowplot"}: a \code{ggplot} grob (via \code{cowplot}).
#'   \item For \code{engine = "patchwork"}: a \code{patchwork} object.
#'   \item For \code{engine = "gridExtra"}: a \code{grid} grob.
#' }
#'
#' @section Dependencies:
#' Requires \pkg{ggplot2}. Depending on \code{engine}, also requires one of
#' \pkg{cowplot}, \pkg{patchwork}, or \pkg{gridExtra}. These are checked at runtime.
#'
#' @examples
#' \dontrun{
#' my_plot_management_grid(fit)
#' my_plot_management_grid(fit, engine = "patchwork", title = "Assessment summary")
#' my_plot_management_grid(fit, CI = 0.8)
#' }
#'
#' @import ggplot2
#' @export
my_plot_management_grid_New <- function(rep,
                                        title = NULL,
                                        widths = c(1, 1, 1),
                                        heights = c(1, 1),
                                        engine = c("cowplot", "patchwork", "gridExtra"),
                                        patchwork_guides = c("keep", "collect"),
                                        ...) {
  # ---- input checks / normalization ----
  if (!inherits(rep, "spictcls")) {
    stop("`rep` must be a fitted SPiCT report object (class 'spictcls').")
  }
  engine <- match.arg(engine)
  patchwork_guides <- match.arg(patchwork_guides)
  dots <- list(...)

  # normalize widths/heights to length 3 and 2 respectively
  widths  <- as.numeric(widths)
  heights <- as.numeric(heights)
  if (!length(widths))  widths <- c(1,1,1)
  if (!length(heights)) heights <- c(1,1)
  if (length(widths)  != 3) widths  <- rep_len(widths, 3)
  if (length(heights) != 2) heights <- rep_len(heights, 2)
  if (!all(is.finite(widths)) || !all(widths > 0))   widths  <- c(1,1,1)
  if (!all(is.finite(heights)) || !all(heights > 0)) heights <- c(1,1)

  title <- if (is.null(title) || !nzchar(as.character(title))) NULL else as.character(title)

  # ---- tiny internals ----
  .fn_get <- function(name) {
    f <- get0(name, mode = "function", inherits = TRUE)
    if (is.null(f)) {
      # Return a dummy function that throws; caller will catch & draw placeholder.
      force(name)
      return(function(...) stop(sprintf("Required function '%s' not found.", name)))
    }
    f
  }

  .filter_args <- function(fun, dots) {
    # Only forward args that the function actually knows (excluding "...")
    fml <- names(formals(fun))
    if (is.null(fml)) return(list())
    fml <- setdiff(fml, "...")
    keep <- intersect(names(dots), fml)
    if (!length(keep)) return(list())
    dots[keep]
  }

  .placeholder_plot <- function(label, err = NULL) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.6, label = label,
                        hjust = 0.5, vjust = 0.5, size = 4, fontface = "bold") +
      ggplot2::annotate("text", x = 0.5, y = 0.4,
                        label = if (is.null(err)) "Panel unavailable." else paste("Error:", err),
                        hjust = 0.5, vjust = 0.5, size = 3) +
      ggplot2::theme_void()
  }

  .safe_make <- function(expr, label) {
    out <- try(expr(), silent = TRUE)
    if (inherits(out, "try-error")) {
      return(.placeholder_plot(label, err = as.character(attr(out, "condition")$message)))
    }
    # Accept ggplot, grob, gtable, patchwork; otherwise coerce or fallback
    if (inherits(out, "ggplot")) return(out)
    if (inherits(out, c("grob", "gTree", "gList", "gtable"))) return(out)
    # Some user functions may return patchwork objects
    if (inherits(out, "patchwork")) return(out)
    # Last resort: draw a placeholder to avoid crashing layout engines
    .placeholder_plot(label, err = "Unrecognized panel return type.")
  }

  # When engines need grobs, coerce ggplot/patchwork to grob safely
  .as_grob <- function(x) {
    if (inherits(x, c("grob", "gTree", "gList", "gtable"))) return(x)
    if (inherits(x, "ggplot")) return(ggplot2::ggplotGrob(x))
    # For patchwork: convert by drawing to grob via grid::convertWidth/Height trick
    if (inherits(x, "patchwork")) {
      if (!requireNamespace("patchwork", quietly = TRUE)) {
        return(.placeholder_plot("Panel", "patchwork not installed; cannot convert."))
      }
      # patchwork prints as ggplot-like; grab a grob from the built object
      p <- patchwork::wrap_plots(x)
      return(ggplot2::ggplotGrob(p))
    }
    # Unknown type
    .placeholder_plot("Panel", "Cannot coerce to grob.")
  }

  # ---- panel factories ----
  f_B    <- .fn_get("my_plot_manage_biomass_panel")
  f_F    <- .fn_get("my_plot_manage_f_panel")
  f_C    <- .fn_get("my_plot_manage_catch_panel")
  f_BB   <- .fn_get("my_plot_manage_bbmsy_panel")
  f_FF   <- .fn_get("my_plot_manage_ffmsy_panel")
  f_KOBE <- .fn_get("my_plot_kobe_all_management_scenario")

  has_man <- ("man" %in% names(rep)) && length(rep$man) > 0

  # Base args per panel (respect caller's ... only where relevant)
  args_B  <- c(list(rep_man = rep), .filter_args(f_B,  dots))
  args_F  <- c(list(rep_man = rep), .filter_args(f_F,  dots))
  args_C  <- c(list(rep_man = rep), .filter_args(f_C,  dots))
  args_BB <- c(list(rep_man = rep), .filter_args(f_BB, dots))
  args_FF <- c(list(rep_man = rep), .filter_args(f_FF, dots))
  args_K  <- c(list(rep = rep),     .filter_args(f_KOBE, dots))

  # Legend policy for cowplot when scenarios exist:
  # show on Biomass only (if that panel happens to accept `show_legend`), hide elsewhere.
  if (engine == "cowplot" && has_man) {
    if (!("show_legend" %in% names(args_B)))  args_B$show_legend  <- TRUE
    if (!("show_legend" %in% names(args_F)))  args_F$show_legend  <- FALSE
    if (!("show_legend" %in% names(args_C)))  args_C$show_legend  <- FALSE
    if (!("show_legend" %in% names(args_BB))) args_BB$show_legend <- FALSE
    if (!("show_legend" %in% names(args_FF))) args_FF$show_legend <- FALSE
  }

  # ---- build panels (fault-tolerant) ----
  p1 <- .safe_make(function() do.call(f_B,    args_B),  "Biomass panel")
  p2 <- .safe_make(function() do.call(f_F,    args_F),  "Fishing mortality")
  p3 <- .safe_make(function() do.call(f_C,    args_C),  "Catch")
  p4 <- .safe_make(function() do.call(f_BB,   args_BB), "B/Bmsy")
  p5 <- .safe_make(function() do.call(f_FF,   args_FF), "F/Fmsy")
  p6 <- .safe_make(function() do.call(f_KOBE, args_K),  "Kobe")

  # ---- layout engines ----
  if (engine == "cowplot") {
    if (!requireNamespace("cowplot", quietly = TRUE)) {
      stop("Please install 'cowplot' for engine = 'cowplot'.")
    }
    # cowplot tolerates both ggplot and grob via cowplot::as_grob
    g1 <- cowplot::as_grob(p1); g2 <- cowplot::as_grob(p2); g3 <- cowplot::as_grob(p3)
    g4 <- cowplot::as_grob(p4); g5 <- cowplot::as_grob(p5); g6 <- cowplot::as_grob(p6)

    row1 <- cowplot::plot_grid(g1, g2, g3, ncol = 3, rel_widths = widths, align = "hv", axis = "tblr")
    row2 <- cowplot::plot_grid(g4, g5, g6, ncol = 3, rel_widths = widths, align = "hv", axis = "tblr")
    out  <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = heights)
    if (!is.null(title)) {
      out <- cowplot::ggdraw() +
        cowplot::draw_plot(out) +
        cowplot::draw_label(title, x = 0.5, y = 0.98, hjust = 0.5, vjust = 1,
                            fontface = "bold", size = 14)
    }
    return(out)
  }

  if (engine == "patchwork") {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Please install 'patchwork' for engine = 'patchwork'.")
    }
    # Patchwork expects ggplot-like objects; coerce grobs to ggplot where possible
    to_plot <- function(x) {
      if (inherits(x, "ggplot")) return(x)
      if (inherits(x, "patchwork")) return(x)
      # Convert grob -> ggplot via ggplotify if available; else wrap in placeholder
      if (inherits(x, c("grob", "gTree", "gList", "gtable"))) {
        if (requireNamespace("ggplotify", quietly = TRUE)) {
          return(ggplotify::as.ggplot(x))
        } else {
          return(.placeholder_plot("Panel", "Install 'ggplotify' to render non-ggplot panels in patchwork."))
        }
      }
      .placeholder_plot("Panel", "Unsupported panel type for patchwork.")
    }

    grid_obj <- (to_plot(p1) | to_plot(p2) | to_plot(p3)) /
      (to_plot(p4) | to_plot(p5) | to_plot(p6))

    grid_obj <- grid_obj + patchwork::plot_layout(
      widths = widths,
      heights = heights,
      guides = if (patchwork_guides == "collect") "collect" else "keep"
    )
    if (!is.null(title)) {
      grid_obj <- grid_obj + patchwork::plot_annotation(
        title = title,
        theme = ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
      )
    }
    return(grid_obj)
  }

  # gridExtra
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Please install 'gridExtra' for engine = 'gridExtra'.")
  }
  # gridExtra works on grobs; coerce everything
  g1 <- .as_grob(p1); g2 <- .as_grob(p2); g3 <- .as_grob(p3)
  g4 <- .as_grob(p4); g5 <- .as_grob(p5); g6 <- .as_grob(p6)

  row1 <- gridExtra::arrangeGrob(g1, g2, g3, ncol = 3, widths = widths)
  row2 <- gridExtra::arrangeGrob(g4, g5, g6, ncol = 3, widths = widths)
  gridExtra::arrangeGrob(row1, row2, ncol = 1, heights = heights,
                         top = if (!is.null(title)) title else NULL)
}
