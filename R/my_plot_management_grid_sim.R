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
#' Absolute Biomass panel by default and suppressed on the other five panels; the
#' Kobe plot keeps its own legend behavior.
#'
#' @param rep A fitted SPiCT report object (\code{spictcls}). May include
#'   management scenarios in \code{rep$man}.
#' @param title Optional title drawn above the whole grid (engine dependent).
#' @param widths Numeric vector of length 3 with relative column widths. Default \code{c(1,1,1)}.
#' @param heights Numeric vector of length 2 with relative row heights. Default \code{c(1,1)}.
#' @param engine Which layout backend to use: \code{"cowplot"}, \code{"patchwork"},
#'   or \code{"gridExtra"}. Default \code{"cowplot"}.
#' @param patchwork_guides For \code{engine = "patchwork"}, whether to \code{"keep"}
#'   per-panel guides or \code{"collect"} them. Default \code{"keep"}.
#' @param ... Additional arguments passed through to the underlying panel functions.
#'   Only arguments matching a panel’s formal parameters are forwarded to that panel.
#'
#' @details
#' Panels used (looked up at runtime; errors render a labeled placeholder panel):
#' \itemize{
#'   \item \code{my_plot_manage_biomass_panel_sim(rep_man = rep, ...)}
#'   \item \code{my_plot_manage_f_panel_sim(rep_man = rep, ...)}
#'   \item \code{my_plot_manage_catch_panel_sim(rep_man = rep, ...)}
#'   \item \code{my_plot_manage_bbmsy_panel_sim(rep_man = rep, ...)}
#'   \item \code{my_plot_manage_ffmsy_panel_sim(rep_man = rep, ...)}
#'   \item \code{my_plot_kobe_all_management_scenario_sim(rep = rep, ...)}
#' }
#'
#' The function keeps each panel’s own settings (titles, themes, guides) unless
#' you override them via \code{...}. For \code{engine = "cowplot"} with scenarios,
#' the function enforces a clear default: scenario legend shown on **Biomass** only,
#' hidden on the other non-Kobe panels; Kobe’s legend remains as defined by
#' \code{my_plot_kobe_all_management_scenario_sim()}.
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
#' # Minimal:
#' my_plot_management_grid_sim(fit)
#'
#' # With patchwork and a title:
#' my_plot_management_grid_sim(fit, engine = "patchwork",
#'                             title = "Assessment summary")
#'
#' # Forward an argument only some panels know (e.g., CI level):
#' my_plot_management_grid_sim(fit, CI = 0.8)
#' }
#'
#' @seealso
#' \code{\link{my_plot_manage_biomass_panel_sim}},
#' \code{\link{my_plot_manage_f_panel_sim}},
#' \code{\link{my_plot_manage_catch_panel_sim}},
#' \code{\link{my_plot_manage_bbmsy_panel_sim}},
#' \code{\link{my_plot_manage_ffmsy_panel_sim}},
#' \code{\link{my_plot_kobe_all_management_scenario_sim}}
#'
#' @import ggplot2
#' @export
my_plot_management_grid_sim <- function(rep,
                                        title = NULL,
                                        widths = c(1, 1, 1),
                                        heights = c(1, 1),
                                        engine = c("cowplot", "patchwork", "gridExtra"),
                                        patchwork_guides = c("keep", "collect"),
                                        ...) {
  stopifnot(inherits(rep, "spictcls"))
  engine <- match.arg(engine)
  patchwork_guides <- match.arg(patchwork_guides)
  dots <- list(...)

  .fn_get <- function(name) {
    f <- get0(name, mode = "function", inherits = TRUE)
    if (is.null(f)) stop(sprintf("Required function '%s' not found in scope.", name))
    f
  }
  .filter_args <- function(fun, dots) {
    fml <- names(formals(fun)); fml <- setdiff(fml, "...")
    keep <- intersect(names(dots), fml); if (!length(keep)) return(list())
    dots[keep]
  }
  .safe_plot <- function(make_plot_expr, label_when_error) {
    tryCatch(make_plot_expr(), error = function(e) {
      ggplot2::ggplot() +
        ggplot2::annotate("text", x=.5, y=.6, label=label_when_error,
                          hjust=.5, vjust=.5, size=4, fontface="bold") +
        ggplot2::annotate("text", x=.5, y=.4, label=paste("Error:", e$message),
                          hjust=.5, vjust=.5, size=3) +
        ggplot2::theme_void()
    })
  }

  # functions
  f_B    <- .fn_get("my_plot_manage_biomass_panel_sim")
  f_F    <- .fn_get("my_plot_manage_f_panel_sim")
  f_C    <- .fn_get("my_plot_manage_catch_panel_sim")
  f_BB   <- .fn_get("my_plot_manage_bbmsy_panel_sim")
  f_FF   <- .fn_get("my_plot_manage_ffmsy_panel_sim")
  f_KOBE <- .fn_get("my_plot_kobe_all_management_scenario_sim")

  # detect scenarios
  has_man <- ("man" %in% names(rep)) && length(rep$man) > 0

  # Build per-plot call args (respect caller's ... when relevant)
  args_B  <- c(list(rep_man = rep), .filter_args(f_B,  dots))
  args_F  <- c(list(rep_man = rep), .filter_args(f_F,  dots))
  args_C  <- c(list(rep_man = rep), .filter_args(f_C,  dots))
  args_BB <- c(list(rep_man = rep), .filter_args(f_BB, dots))
  args_FF <- c(list(rep_man = rep), .filter_args(f_FF, dots))
  args_K  <- c(list(rep = rep),     .filter_args(f_KOBE, dots))

  # --- Legend fix for cowplot with scenarios present ---
  # Show scenario legend ONLY on Absolute biomass panel.
  if (engine == "cowplot" && has_man) {
    # Only set if not explicitly provided by the caller
    if (!"show_legend" %in% names(args_B))  args_B$show_legend  <- TRUE
    if (!"show_legend" %in% names(args_F))  args_F$show_legend  <- FALSE
    if (!"show_legend" %in% names(args_C))  args_C$show_legend  <- FALSE
    if (!"show_legend" %in% names(args_BB)) args_BB$show_legend <- FALSE
    if (!"show_legend" %in% names(args_FF)) args_FF$show_legend <- FALSE
    # Kobe stays "as is" (your function already controls its legend internally)
  }

  # make plots
  p1 <- .safe_plot(function() do.call(f_B,    args_B),  "Biomass panel")
  p2 <- .safe_plot(function() do.call(f_F,    args_F),  "Fishing mortality")
  p3 <- .safe_plot(function() do.call(f_C,    args_C),  "Catch")
  p4 <- .safe_plot(function() do.call(f_BB,   args_BB), "B/Bmsy")
  p5 <- .safe_plot(function() do.call(f_FF,   args_FF), "F/Fmsy")
  p6 <- .safe_plot(function() do.call(f_KOBE, args_K),  "Kobe")

  if (engine == "cowplot") {
    if (!requireNamespace("cowplot", quietly = TRUE))
      stop("Please install 'cowplot' for engine='cowplot'.")
    row1 <- cowplot::plot_grid(p1, p2, p3, ncol = 3, rel_widths = widths, align = "hv", axis = "tblr")
    row2 <- cowplot::plot_grid(p4, p5, p6, ncol = 3, rel_widths = widths, align = "hv", axis = "tblr")
    out  <- cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = heights)
    if (!is.null(title) && nzchar(title)) {
      out <- cowplot::ggdraw() +
        cowplot::draw_plot(out) +
        cowplot::draw_label(title, x = 0.5, y = 0.98, hjust = 0.5, vjust = 1,
                            fontface = "bold", size = 14)
    }
    return(out)
  }

  if (engine == "patchwork") {
    if (!requireNamespace("patchwork", quietly = TRUE))
      stop("Please install 'patchwork' for engine='patchwork'.")
    grid_obj <- (p1 | p2 | p3) / (p4 | p5 | p6)
    grid_obj <- grid_obj + patchwork::plot_layout(widths = widths, heights = heights,
                                                  guides = if (patchwork_guides == "collect") "collect" else "keep")
    if (!is.null(title) && nzchar(title)) {
      grid_obj <- grid_obj + patchwork::plot_annotation(title = title,
                                                        theme = ggplot2::theme(plot.margin = ggplot2::margin(0,0,0,0)))
    }
    return(grid_obj)
  }

  # gridExtra fallback
  if (!requireNamespace("gridExtra", quietly = TRUE))
    stop("Please install 'gridExtra' for engine='gridExtra'.")
  row1 <- gridExtra::arrangeGrob(p1, p2, p3, ncol = 3, widths = widths)
  row2 <- gridExtra::arrangeGrob(p4, p5, p6, ncol = 3, widths = widths)
  grob <- gridExtra::arrangeGrob(row1, row2, ncol = 1, heights = heights,
                                 top = if (!is.null(title) && nzchar(title)) title else NULL)
  grob
}
