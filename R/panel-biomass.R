# ===========================
# 1/6 — Biomass panel
# ===========================

#' Biomass panel (with Bmsy band and B/Bmsy-on-B ribbon)
#'
#' Plots estimated \eqn{B_t} with dotted CIs, overlays the \eqn{B_{MSY}}
#' band/line and a blue \eqn{B/B_{MSY}} ribbon on the \eqn{B} scale, adds a
#' vertical line at end-of-observations, and optionally shows q-scaled index points.
#' A secondary y-axis displays \eqn{B_t/B_{MSY}}.
#'
#' @param model A fitted SPiCT object (\code{spictcls}).
#' @param line_color Colour for the estimate line.
#' @param show_CIs Logical; draw ribbons and dotted CI edges.
#' @param CI Confidence level for intervals (default \code{0.95}).
#'
#' @return A \code{ggplot} object.
#' @family ELU2-panels
#' @seealso \code{\link{plot_elu2_panel_bbmsy}}, \code{\link{plot_elu2_panel_f}},
#'   \code{\link{plot_elu2_panel_ffmsy}}, \code{\link{plot_elu2_panel_catch}},
#'   \code{\link{plot_elu2_panel_kobe}}
#' @export
#' @examples
#' \dontrun{
#' m1 <- all_models$S1$S1P.SDM
#' plot_elu2_panel_biomass(m1)
#' }
plot_elu2_panel_biomass <- function(model,
                                    line_color = "blue",
                                    show_CIs   = TRUE,
                                    CI         = 0.95) {

  stopifnot(inherits(model, "spictcls"))
  vline_col  <- "grey50"; vline_size <- 0.2
  manflag <- ("man" %in% names(model))

  Best <- get.par("logB",     model, exp = TRUE, CI = CI)
  BB   <- get.par("logBBmsy", model, exp = TRUE, CI = CI)

  repmax   <- if (manflag) get.manmax(model) else model
  Bmsy_all <- get.par("logBmsy", repmax, exp = TRUE, CI = CI)
  Bmsyvec  <- get.msyvec(repmax$inp, Bmsy_all)
  Bmsy     <- if (!is.null(nrow(Bmsy_all))) Bmsy_all[1, ] else Bmsy_all

  df_B <- data.frame(
    time = .spict_time_from_par(model, Best),
    lwr  = Best[, 1], est = Best[, 2], upr = Best[, 3]
  )

  inp <- model$inp
  ind_in <- inp$indest
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(length(ind_in) - 1)]
  ind_pr <- if (!manflag) inp$indpred else integer(0)

  df_B_in <- if (length(ind_in)) df_B[ind_in, , drop = FALSE] else df_B[0, ]
  df_B_pr <- if (length(ind_pr)) df_B[ind_pr, , drop = FALSE] else df_B[0, ]

  Bmsy_est <- Bmsy[2]
  df_BB_rib <- .spict_make_ribbon_df(
    time = .spict_time_from_par(model, BB),
    lwr  = BB[, 1] * Bmsy_est,
    upr  = BB[, 3] * Bmsy_est
  )

  df_Bmsy_band <- data.frame(time = repmax$inp$time, ymin = Bmsyvec$ll, ymax = Bmsyvec$ul)
  df_Bmsy_line <- data.frame(time = repmax$inp$time, y = Bmsyvec$msy)

  p <- ggplot2::ggplot()
  if (isTRUE(show_CIs)) {
    p <- p + ggplot2::geom_ribbon(data = df_Bmsy_band, ggplot2::aes(x = time, ymin = ymin, ymax = ymax), fill = "grey80")
    if (nrow(df_BB_rib)) {
      p <- p + ggplot2::geom_ribbon(data = df_BB_rib, ggplot2::aes(x = time, ymin = ymin, ymax = ymax), fill = grDevices::rgb(0,0,1,0.10)) +
        ggplot2::geom_line(data = transform(df_BB_rib, y = ymin), ggplot2::aes(x = time, y = y), color = grDevices::rgb(0,0,1,0.20), linewidth = 0.6) +
        ggplot2::geom_line(data = transform(df_BB_rib, y = ymax), ggplot2::aes(x = time, y = y), color = grDevices::rgb(0,0,1,0.20), linewidth = 0.6)
    }
  }
  p <- p + ggplot2::geom_line(data = df_Bmsy_line, ggplot2::aes(x = time, y = y), color = "black", linewidth = 0.7)

  if (isTRUE(show_CIs) && nrow(df_B_in)) {
    p <- p + ggplot2::geom_line(data = df_B_in, ggplot2::aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
      ggplot2::geom_line(data = df_B_in, ggplot2::aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
  }
  if (nrow(df_B_in)) p <- p + ggplot2::geom_line(data = df_B_in, ggplot2::aes(x = time, y = est), linewidth = 0.8, color = line_color)
  if (!manflag && nrow(df_B_pr)) {
    p <- p + ggplot2::geom_line(data = df_B_pr, ggplot2::aes(x = time, y = est), linetype = "dashed", linewidth = 0.8, color = line_color)
    if (isTRUE(show_CIs)) {
      p <- p + ggplot2::geom_line(data = df_B_pr, ggplot2::aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
        ggplot2::geom_line(data = df_B_pr, ggplot2::aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
    }
  }

  obs_end <- .spict_obs_end_overall(model)
  if (is.finite(obs_end)) p <- p + ggplot2::geom_vline(xintercept = obs_end, color = vline_col, linewidth = vline_size)

  # q-scaled index points
  idx <- .spict_index_points(model)
  if (length(idx) >= 1) p <- p + ggplot2::geom_point(data = idx[[1]], ggplot2::aes(x = time, y = obs), color = "blue", shape = 16, size = 2, inherit.aes = FALSE)
  if (length(idx) >= 2) p <- p + ggplot2::geom_point(data = idx[[2]], ggplot2::aes(x = time, y = obs), shape = 22, color = "black", fill = "green", size = 2, stroke = 0.5, inherit.aes = FALSE)

  p <- p +
    ggplot2::labs(title = "Absolute biomass", x = "Year", y = expression(B[t])) +
    .spict_theme_minimal_compact2() +
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ . / Bmsy_est, name = expression(B[t]/B[MSY])))

  p

}  # (body unchanged). Update and give the final version
