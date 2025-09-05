# ===========================
# 5/6 — F/Fmsy panel
# ===========================

#' F/Fmsy panel (ribbon + outlines; seasonal overlay if applicable)
#'
#' Plots \eqn{F_t/F_{MSY}} with a light blue ribbon and dotted CI outlines,
#' solid line over data, dashed over predictions, reference line at 1, and
#' a vertical line at end-of-observations. If seasonal, overlays sub-annual
#' \eqn{F_t/F_{MSY}}.
#'
#' @param model A fitted SPiCT object (\code{spictcls}).
#' @param line_color Colour for the estimate line.
#' @param show_CIs Logical; draw ribbons and dotted CI edges.
#' @param CI Confidence level for intervals (default \code{0.95}).
#'
#' @return A \code{ggplot} object.
#' @family ELU2-panels
#' @export
#' @examples
#' \dontrun{
#' m1 <- all_models$S1$S1P.SDM
#' plot_elu2_panel_ffmsy(m1)
#' }
plot_elu2_panel_ffmsy <- function(model,
                                  line_color = "#1b9e77",
                                  show_CIs   = TRUE,
                                  CI         = 0.95) {
  stopifnot(inherits(model, "spictcls"))
  vline_col  <- "grey50"; vline_size <- 0.2
  manflag <- ("man" %in% names(model))
  inp <- model$inp

  FF <- get.par("logFFmsy", model, exp = TRUE, CI = CI)
  df <- data.frame(time = .spict_time_from_par(model, FF), lwr = FF[,1], est = FF[,2], upr = FF[,3])

  ind_in <- inp$indest
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(length(ind_in) - 1)]
  ind_pr <- if (!manflag) inp$indpred else integer(0)

  df_in <- if (length(ind_in)) df[ind_in, , drop = FALSE] else df[0, ]
  df_pr <- if (length(ind_pr)) df[ind_pr, , drop = FALSE] else df[0, ]

  rib <- .spict_make_ribbon_df(df$time, df$lwr, df$upr)

  p <- ggplot2::ggplot()
  if (isTRUE(show_CIs) && nrow(rib)) {
    p <- p + ggplot2::geom_ribbon(data = rib, ggplot2::aes(x = time, ymin = ymin, ymax = ymax), fill = grDevices::rgb(0,0,1,0.10)) +
      ggplot2::geom_line(data = transform(rib, y = ymin), ggplot2::aes(x = time, y = y), color = grDevices::rgb(0,0,1,0.20), linewidth = 0.6) +
      ggplot2::geom_line(data = transform(rib, y = ymax), ggplot2::aes(x = time, y = y), color = grDevices::rgb(0,0,1,0.20), linewidth = 0.6)
  }
  if (nrow(df_in)) p <- p + ggplot2::geom_line(data = df_in, ggplot2::aes(x = time, y = est), color = line_color, linewidth = 0.8)
  if (nrow(df_pr)) p <- p + ggplot2::geom_line(data = df_pr, ggplot2::aes(x = time, y = est), linetype = "dashed", color = line_color, linewidth = 0.8)

  obs_end <- .spict_obs_end_overall(model)
  if (is.finite(obs_end)) p <- p + ggplot2::geom_vline(xintercept = obs_end, color = vline_col, linewidth = vline_size)

  p <- p +
    ggplot2::geom_hline(yintercept = 1, linetype = "solid", color = "black", linewidth = 0.8) +
    ggplot2::labs(title = "Relative fishing mortality", x = "Year", y = expression(bold(F/F[MSY]))) +
    .spict_theme_minimal_compact2()

  # Seasonal overlay (unchanged)
  if (min(inp$dtc, na.rm = TRUE) < 1 && nrow(df_in)) {
    FFs <- get.par("logFFmsy", model, exp = TRUE, CI = CI)
    p <- p + ggplot2::geom_line(
      data = data.frame(time = df_in$time, est = FFs[seq_len(nrow(df_in)), 2]),
      ggplot2::aes(x = time, y = est),
      color = grDevices::rgb(0,0,1,0.4), linewidth = 0.6
    )
  }

  p

}  # (body unchanged)
