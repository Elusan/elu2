# ===========================
# 4/6 — Absolute Ft panel
# ===========================

#' Absolute Ft panel (median + CI; Fmsy band + line; optional sub-annual Ft)
#'
#' Plots absolute \eqn{F_t} with dotted CI outlines, overlays the \eqn{F_{MSY}}
#' band/line, optionally shows sub-annual \eqn{F_t} (seasonal models),
#' and adds a vertical line at end-of-observations. Secondary y-axis can show
#' \eqn{F_t/F_{MSY}} when \eqn{F_{MSY}} is constant.
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
#' plot_elu2_panel_f(m1)
#' }
plot_elu2_panel_f <- function(model,
                              line_color = "#1b9e77",
                              show_CIs   = TRUE,
                              CI         = 0.95) {

  stopifnot(inherits(model, "spictcls"))
  vline_col  <- "grey50"; vline_size <- 0.2
  manflag <- ("man" %in% names(model))
  inp <- model$inp

  repmax  <- if (manflag) get.manmax(model) else model
  tvgflag <- isTRUE(repmax$inp$timevaryinggrowth) || isTRUE(repmax$inp$logmcovflag)

  Fest <- get.par("logFnotS", model, exp = TRUE, CI = CI)
  dfF <- data.frame(
    time = .spict_time_from_par(model, Fest),
    lwr  = Fest[, 1], est = Fest[, 2], upr = Fest[, 3]
  )

  ind_in <- inp$indest
  if (length(ind_in) > 0) ind_in <- ind_in[seq_len(length(ind_in) - 1)]
  ind_pr <- if (!manflag) inp$indpred else integer(0)
  dfF_in <- if (length(ind_in)) dfF[ind_in, , drop = FALSE] else dfF[0, ]
  dfF_pr <- if (length(ind_pr)) dfF[ind_pr, , drop = FALSE] else dfF[0, ]

  if (tvgflag) {
    Fmsy_all <- get.par("logFmsyvec", repmax, exp = TRUE, CI = CI)
    Fmsyvec  <- as.data.frame(Fmsy_all); Fmsyvec$msy <- Fmsyvec$est
    Fmsy_band <- data.frame(time = repmax$inp$time, ymin = Fmsyvec$ll, ymax = Fmsyvec$ul)
    Fmsy_line <- data.frame(time = repmax$inp$time, y = Fmsyvec$msy)
    fmsy_mult <- matrix(rep(Fmsyvec$msy, each = 3), ncol = 3, byrow = TRUE)
    sec_axis_obj <- ggplot2::waiver()
  } else {
    Fmsy_all <- get.par("logFmsy", repmax, exp = TRUE, CI = CI)
    if (any(is.na(Fmsy_all))) Fmsy_all <- get.par("logFmsyd", repmax, exp = TRUE, CI = CI)
    xm <- as.matrix(Fmsy_all)
    col_est <- if (!is.null(colnames(xm))) which(colnames(xm) %in% c("est","Est","EST"))[1] else NA_integer_
    if (is.na(col_est)) col_est <- 2
    Fmsy_est <- as.numeric(xm[1, col_est])

    Fmsyvec  <- get.msyvec(repmax$inp, Fmsy_all)
    Fmsy_band <- data.frame(time = repmax$inp$time, ymin = Fmsyvec$ll, ymax = Fmsyvec$ul)
    Fmsy_line <- data.frame(time = repmax$inp$time, y = Fmsyvec$msy)
    fmsy_mult <- Fmsy_est
    sec_axis_obj <- ggplot2::sec_axis(~ . / Fmsy_est, name = expression(F[t]/F[MSY]))
  }

  FFrel <- get.par("logFFmsynotS", model, exp = TRUE, CI = CI)[, 1:3]
  last_obs_time <- if (manflag) repmax$inp$timerange[2] else max(repmax$inp$time, na.rm = TRUE)
  ix_last <- which(repmax$inp$time == last_obs_time)
  if (!length(ix_last)) ix_last <- length(repmax$inp$time)
  indxmax <- max(ix_last - if (manflag) 1 else 0, 1)

  timef <- repmax$inp$time[seq_len(indxmax)]
  if (is.matrix(fmsy_mult)) {
    FFabs <- FFrel * fmsy_mult
  } else {
    FFabs <- FFrel * as.numeric(fmsy_mult)
  }
  ribFF <- .spict_make_ribbon_df(timef, FFabs[seq_len(indxmax), 1], FFabs[seq_len(indxmax), 3])

  p <- ggplot2::ggplot()
  if (isTRUE(show_CIs)) {
    p <- p + ggplot2::geom_ribbon(data = Fmsy_band, ggplot2::aes(x = time, ymin = ymin, ymax = ymax), fill = "grey80")
  }
  p <- p + ggplot2::geom_line(data = Fmsy_line, ggplot2::aes(x = time, y = y), color = "black", linewidth = 0.7)

  if (isTRUE(show_CIs) && nrow(ribFF)) {
    p <- p + ggplot2::geom_ribbon(data = ribFF, ggplot2::aes(x = time, ymin = ymin, ymax = ymax), fill = grDevices::rgb(0,0,1,0.10)) +
      ggplot2::geom_line(data = transform(ribFF, y = ymin), ggplot2::aes(x = time, y = y), color = grDevices::rgb(0,0,1,0.20), linewidth = 0.6) +
      ggplot2::geom_line(data = transform(ribFF, y = ymax), ggplot2::aes(x = time, y = y), color = grDevices::rgb(0,0,1,0.20), linewidth = 0.6)
  }

  if (min(inp$dtc, na.rm = TRUE) < 1) {
    sF <- get.par("logFs", model, exp = TRUE, CI = CI)
    p <- p + ggplot2::geom_line(
      data = data.frame(time = repmax$inp$time[seq_len(indxmax)], est = sF[seq_len(indxmax), 2]),
      ggplot2::aes(x = time, y = est),
      color = grDevices::rgb(0,0,1,0.4), linewidth = 0.6
    )
  }

  if (isTRUE(show_CIs) && nrow(dfF_in)) {
    p <- p + ggplot2::geom_line(data = dfF_in, ggplot2::aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
      ggplot2::geom_line(data = dfF_in, ggplot2::aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
  }
  if (nrow(dfF_in)) p <- p + ggplot2::geom_line(data = dfF_in, ggplot2::aes(x = time, y = est), color = line_color, linewidth = 0.8)

  if (!manflag && nrow(dfF_pr)) {
    if (isTRUE(show_CIs)) {
      p <- p + ggplot2::geom_line(data = dfF_pr, ggplot2::aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
        ggplot2::geom_line(data = dfF_pr, ggplot2::aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
    }
    p <- p + ggplot2::geom_line(data = dfF_pr, ggplot2::aes(x = time, y = est), linetype = "dashed", color = line_color, linewidth = 0.8)
  }

  if (!is.null(inp$timeE) && length(inp$timeE) && !is.null(inp$obsE) && !is.null(inp$dte)) {
    qf <- try(get.par("logqf", model, exp = TRUE, CI = CI), silent = TRUE)
    if (!inherits(qf, "try-error") && is.finite(qf[2])) {
      df_obsF <- data.frame(time = inp$timeE, val = inp$obsE / inp$dte * qf[2])
      p <- p + ggplot2::geom_point(data = df_obsF, ggplot2::aes(x = time, y = val), color = "black", size = 1.8)
    }
  }

  obs_end <- .spict_obs_end_overall(model)
  if (is.finite(obs_end)) p <- p + ggplot2::geom_vline(xintercept = obs_end, color = vline_col, linewidth = vline_size)

  p <- p + ggplot2::labs(title = "Absolute fishing mortality", x = "Year", y = expression(F[t])) +
    .spict_theme_minimal_compact2() +
    ggplot2::scale_y_continuous(sec.axis = sec_axis_obj)

  p

}  # (body unchanged)
