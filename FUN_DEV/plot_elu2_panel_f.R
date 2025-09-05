#' Absolute Ft panel (median + CI; Fmsy band + line; optional sub-annual Ft)
#'
#' @param model A fitted SPiCT object (`spictcls`).
#' @param line_color Line color for estimates.
#' @param show_CIs Logical; draw ribbons & dotted CI edges.
#' @param CI Confidence level for intervals (default 0.95).
#' @return A ggplot object.
#' @export
#' @import ggplot2
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

  p <- ggplot()
  if (isTRUE(show_CIs)) {
    p <- p + geom_ribbon(data = Fmsy_band, aes(x = time, ymin = ymin, ymax = ymax), fill = "grey80")
  }
  p <- p + geom_line(data = Fmsy_line, aes(x = time, y = y), color = "black", linewidth = 0.7)

  if (isTRUE(show_CIs) && nrow(ribFF)) {
    p <- p + geom_ribbon(data = ribFF, aes(x = time, ymin = ymin, ymax = ymax), fill = grDevices::rgb(0,0,1,0.10)) +
      geom_line(data = transform(ribFF, y = ymin), aes(x = time, y = y), color = grDevices::rgb(0,0,1,0.20), linewidth = 0.6) +
      geom_line(data = transform(ribFF, y = ymax), aes(x = time, y = y), color = grDevices::rgb(0,0,1,0.20), linewidth = 0.6)
  }

  if (min(inp$dtc, na.rm = TRUE) < 1) {
    sF <- get.par("logFs", model, exp = TRUE, CI = CI)
    p <- p + geom_line(
      data = data.frame(time = repmax$inp$time[seq_len(indxmax)], est = sF[seq_len(indxmax), 2]),
      aes(x = time, y = est),
      color = grDevices::rgb(0,0,1,0.4), linewidth = 0.6
    )
  }

  if (isTRUE(show_CIs) && nrow(dfF_in)) {
    p <- p + geom_line(data = dfF_in, aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
      geom_line(data = dfF_in, aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
  }
  if (nrow(dfF_in)) p <- p + geom_line(data = dfF_in, aes(x = time, y = est), color = line_color, linewidth = 0.8)

  if (!manflag && nrow(dfF_pr)) {
    if (isTRUE(show_CIs)) {
      p <- p + geom_line(data = dfF_pr, aes(x = time, y = lwr), linetype = "dotted", linewidth = 0.6, color = line_color) +
        geom_line(data = dfF_pr, aes(x = time, y = upr), linetype = "dotted", linewidth = 0.6, color = line_color)
    }
    p <- p + geom_line(data = dfF_pr, aes(x = time, y = est), linetype = "dashed", color = line_color, linewidth = 0.8)
  }

  if (!is.null(inp$timeE) && length(inp$timeE) && !is.null(inp$obsE) && !is.null(inp$dte)) {
    qf <- try(get.par("logqf", model, exp = TRUE, CI = CI), silent = TRUE)
    if (!inherits(qf, "try-error") && is.finite(qf[2])) {
      df_obsF <- data.frame(time = inp$timeE, val = inp$obsE / inp$dte * qf[2])
      p <- p + geom_point(data = df_obsF, aes(x = time, y = val), color = "black", size = 1.8)
    }
  }

  obs_end <- .spict_obs_end_overall(model)
  if (is.finite(obs_end)) p <- p + geom_vline(xintercept = obs_end, color = vline_col, linewidth = vline_size)

  p <- p + labs(title = "Absolute fishing mortality", x = "Year", y = expression(F[t])) +
    .spict_theme_minimal_compact2() +
    scale_y_continuous(sec.axis = sec_axis_obj)

  p
}

# ---- Inlined helpers (only those used by this panel) ----

#' Minimal compact ggplot2 theme (ELU2)
#' @keywords internal
#' @noRd
.spict_theme_minimal_compact2 <- function(base_size = 10, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      axis.text  = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey35", linewidth = 2),
      axis.ticks = ggplot2::element_line(linewidth = 0.5, color = "grey35"),
      axis.ticks.length = grid::unit(3, "pt"),
      strip.background = ggplot2::element_rect(fill = "grey35", color = "grey35", linewidth = 0.5),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      text = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}

#' End of observation time across components
#' @keywords internal
#' @noRd
.spict_obs_end_overall <- function(model) {
  inp <- model$inp
  tr <- inp$timerangeObs
  if (!is.null(tr) && length(tr) >= 2 && is.finite(tr[2])) return(tr[2])
  if (!is.null(inp$timeC) && length(inp$timeC))       return(tail(inp$timeC, 1))
  if (!is.null(inp$timeI) && length(inp$timeI) && length(inp$timeI[[1]]) > 0) return(tail(inp$timeI[[1]], 1))
  if (!is.null(inp$time)  && length(inp$time))        return(max(inp$time, na.rm = TRUE))
  NA_real_
}

#' Build a ribbon dataframe from lwr/upr vectors
#' @keywords internal
#' @noRd
.spict_make_ribbon_df <- function(time, lwr, upr) {
  ok <- is.finite(time) & is.finite(lwr) & is.finite(upr)
  data.frame(time = time[ok], ymin = lwr[ok], ymax = upr[ok])
}

#' Prefer rownames(time) else fallback to inp$time
#' @keywords internal
#' @noRd
.spict_time_from_par <- function(model, par_matrix) {
  rn <- suppressWarnings(as.numeric(rownames(par_matrix)))
  if (length(rn) == nrow(par_matrix) && !all(is.na(rn))) return(rn)
  t_full <- as.numeric(model$inp$time)
  if (length(t_full) >= nrow(par_matrix)) return(t_full[seq_len(nrow(par_matrix))])
  return(rep_len(t_full, nrow(par_matrix)))
}
