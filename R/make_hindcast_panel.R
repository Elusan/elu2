# R/make_hindcast_panel.R

#' Make one hindcast panel (exact ggplot or placeholder)
#'
#' Runs a hindcast on a fitted model and returns a labeled panel using
#' [plotspict.hindcast_elu2_gg_exact_FOR_GRIDS()]. On error, returns a
#' placeholder panel carrying the error note. The label is nudged slightly
#' inward from the top-left border to avoid overlap with the panel frame.
#'
#' @param fit Fitted ELU/SPiCT-style model (no hindcast attached yet).
#' @param panel_label Character to draw at top-left inside the panel (e.g., `"S1P.SDM"`).
#' @param npeels,peel.dtc,mc.cores Hindcast controls; forwarded to [hindcast()].
#' @param CI,verbose,add.mase Passed to [plotspict.hindcast_elu2_gg_exact_FOR_GRIDS()].
#' @param show_legend Logical; keep legend (`TRUE`) or hide (`FALSE`).
#'
#' @return A [ggplot2::ggplot] object.
#' @seealso [hindcast()], [plotspict.hindcast_elu2_gg_exact_FOR_GRIDS()]
#' @export
#' @import ggplot2
make_hindcast_panel <- function(
    fit, panel_label,
    npeels = 2, peel.dtc = FALSE, mc.cores = 1,
    CI = 0.95, verbose = FALSE, add.mase = TRUE,
    show_legend = TRUE
) {
  # tiny inward offsets from the panel border
  add_label <- function(p, lab, h_in = -0.25, v_in = 1.25, lab_size = 3.2) {
    p +
      ggplot2::annotate(
        "text",
        x = -Inf, y = Inf,
        hjust = h_in,  # more negative => a bit more to the right
        vjust = v_in,  # >1 => a bit lower
        label = lab, fontface = 2, size = lab_size
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      (if (!show_legend) ggplot2::theme(legend.position = "none") else ggplot2::theme())
  }

  placeholder <- function(msg = "plot failed; showing placeholder panel.") {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = msg, size = 3) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(6, 6, 6, 8),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey40", linewidth = 1)
      ) |>
      add_label(panel_label)
  }

  out <- try({
    # If you want to be explicit (and avoid clashes), call elu2::hindcast():
    # rep_hc <- elu2::hindcast(fit, npeels = npeels, peel.dtc = peel.dtc, mc.cores = mc.cores)
    rep_hc <- hindcast(fit, npeels = npeels, peel.dtc = peel.dtc, mc.cores = mc.cores)
    p <- plotspict.hindcast_elu2_gg_exact_FOR_GRIDS(
      rep_hc, add.mase = add.mase, CI = CI, verbose = verbose
    )
    add_label(p, panel_label)
  }, silent = TRUE)

  if (inherits(out, "try-error")) placeholder(as.character(out)) else out
}
