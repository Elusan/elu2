# R/theme_minimal_compact2_diag.R
#' Compact diagnostic theme (no frills, bold panels, no legends)
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @return A ggplot2 theme.
#' @export
#' @import ggplot2
theme_minimal_compact2_diag <- function(base_size = 10, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      axis.text  = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "none",
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text  = ggplot2::element_text(size = 14),
      legend.key.size = grid::unit(0.8, "lines"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey45", linewidth = 2),
      axis.ticks = ggplot2::element_line(linewidth = 0.5, color = "grey45"),
      axis.ticks.length = grid::unit(3, "pt"),
      strip.background = ggplot2::element_rect(fill = "grey45", color = "grey45", linewidth = 0.5),
      strip.text = ggplot2::element_text(face = ggplot2::rel(1), size = ggplot2::rel(1)),
      text = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(1, 1, 1, 1)
    )
}

# R/plotspict.diagnostic.process_gg.R
#' Process diagnostics (ggplot2; blue points + light-blue lines on top 4 panels)
#'
#' @export
#' @import ggplot2
#' @importFrom stats acf qnorm ppoints sd
plotspict.diagnostic.process_gg <- function(rep,
                                            lag.max = 4,
                                            qlegend = TRUE,   # kept for API parity; ignored
                                            plot.data = TRUE,
                                            mfcol = FALSE,
                                            add.loess = FALSE,
                                            span = 0.75,
                                            stamp = get.version()) {
  check.rep(rep)
  inp <- rep$inp
  if (!any(names(rep) == "process.resid")) stop("Run calc.process.resid(rep) first.")
  if (!all(c("time","B","F") %in% names(rep$process.resid)))
    stop("rep$process.resid must contain 'time','B','F'.")

  # ---- Aesthetics (BLUE points + LIGHT-BLUE lines) --------------------------
  col_pts  <- "#1f78b4"   # points
  col_line <- "#9ecae1"   # connecting lines
  lwd_main <- 0.5
  pt_size  <- 2.2
  pt_shape <- 19
  zero_lty <- 3

  five_breaks <- function(xmin, xmax) {
    if (!is.finite(xmin) || !is.finite(xmax) || xmin == xmax) return(xmin)
    seq(xmin, xmax, length.out = 5)
  }

  # ---- Diagnostic tests (identical to base via res.diagn) -------------------
  testsB <- res.diagn(rep$process.resid$B, "B", "B")
  testsF <- res.diagn(rep$process.resid$F, "F", "F")
  biasB.p    <- round(testsB$biasB.p,    4)
  biasF.p    <- round(testsF$biasF.p,    4)
  LBoxB.p    <- round(testsB$LBoxB.p,    4)
  LBoxF.p    <- round(testsF$LBoxF.p,    4)
  shapiroB.p <- round(testsB$shapiroB.p, 4)
  shapiroF.p <- round(testsF$shapiroF.p, 4)

  # Title colors matching SPiCT base logic
  colmainB <- ifelse(biasB.p    < 0.05, "red", "forestgreen")  # residuals
  colmainF <- ifelse(biasF.p    < 0.05, "red", "forestgreen")
  colLBoxB <- ifelse(LBoxB.p    < 0.05, "red", "forestgreen")  # ACF
  colLBoxF <- ifelse(LBoxF.p    < 0.05, "red", "forestgreen")
  colShapB <- ifelse(shapiroB.p < 0.05, "red", "forestgreen")  # QQ
  colShapF <- ifelse(shapiroF.p < 0.05, "red", "forestgreen")

  # ---- Time aggregation (same as SPiCT base) --------------------------------
  time.full <- inp$time
  dt <- diff(rep$process.resid$time)[1]
  time.agg <- time.full[which(time.full %% (dt) == 0)]
  time.rep <- rep(time.agg, each = 1 / (inp$dteuler / dt))
  nta <- length(time.agg)

  logB  <- get.par("logB",  rep)[, 2]
  logF  <- get.par("logF",  rep)[, 2]
  logFs <- get.par("logFs", rep)[, 2]
  P     <- rep$obj$report()$P

  logB.agg  <- logF.agg  <- rep(0, nta)
  logFs.agg <- P.agg     <- rep(0, nta)
  for (i in seq_len(nta)) {
    inds <- which(time.agg[i] == time.rep)
    logB.agg[i]  <- mean(logB[inds])
    logF.agg[i]  <- mean(logF[inds])
    logFs.agg[i] <- mean(logFs[inds])
    P.agg[i]     <- sum(P[inds])
  }

  # ---- Season factor (not displayed) ----------------------------------------
  nseas <- inp$nseasons
  q_fac <- floor(((rep$process.resid$time %% 1) * nseas)) + 1L
  q_fac[q_fac < 1L] <- 1L; q_fac[q_fac > nseas] <- nseas
  q_fac <- factor(q_fac, levels = seq_len(nseas))

  # ---- Data frames -----------------------------------------------------------
  d_logB <- data.frame(time = time.agg, logB = logB.agg)
  d_logF <- data.frame(time = time.agg, logF = logF.agg)
  d_resB <- data.frame(time = rep$process.resid$time, res = rep$process.resid$B, seas = q_fac)
  d_resF <- data.frame(time = rep$process.resid$time, res = rep$process.resid$F, seas = q_fac)
  d_naB  <- d_resB[is.na(d_resB$res), , drop = FALSE]
  d_naF  <- d_resF[is.na(d_resF$res), , drop = FALSE]

  # ---- ACF / QQ helpers ------------------------------------------------------
  acf_df <- function(x, lag.max) {
    x <- stats::na.omit(x)
    out <- stats::acf(x, lag.max = lag.max, plot = FALSE, na.action = na.omit)
    data.frame(lag = out$lag[, , 1], acf = out$acf[, , 1], n = length(x))
  }
  d_acfB <- acf_df(rep$process.resid$B, lag.max)
  d_acfF <- acf_df(rep$process.resid$F, lag.max)

  qq_df <- function(x) {
    x <- stats::na.omit(x)
    n <- length(x)
    z <- stats::qnorm(stats::ppoints(n))
    xs <- sort(x)
    mu <- mean(x); sdv <- stats::sd(x)
    data.frame(theoretical = z, sample = xs, ref_x = z, ref_y = mu + sdv * z)
  }
  d_qqB <- qq_df(rep$process.resid$B)
  d_qqF <- qq_df(rep$process.resid$F)

  # ---- Theme (no legends) ----------------------------------------------------
  thm <- theme_minimal_compact2_diag()

  xmin <- min(time.agg, na.rm = TRUE); xmax <- max(time.agg, na.rm = TRUE)
  brk5 <- five_breaks(xmin, xmax)

  # ===== Row 1: logB / logF ===================================================
  p_logB <- ggplot2::ggplot(d_logB, aes(time, logB)) +
    ggplot2::geom_line(linewidth = lwd_main, color = col_line, na.rm = TRUE) +
    ggplot2::geom_point(size = pt_size, shape = pt_shape, color = col_pts, na.rm = TRUE) +
    ggplot2::labs(x = "Time", y = "logB", title = "Biomass") +
    ggplot2::scale_x_continuous(breaks = brk5) +
    thm

  p_logF <- ggplot2::ggplot(d_logF, aes(time, logF)) +
    ggplot2::geom_line(linewidth = lwd_main, color = col_line, na.rm = TRUE) +
    ggplot2::geom_point(size = pt_size, shape = pt_shape, color = col_pts, na.rm = TRUE) +
    ggplot2::labs(x = "Time", y = "logF", title = "Fishing mortality") +
    ggplot2::scale_x_continuous(breaks = brk5) +
    thm

  # ===== Row 2: B/F residuals =================================================
  mk_res_panel <- function(df, df_na, ttl, ttl_col) {
    g <- ggplot2::ggplot(df, aes(time, res)) +
      ggplot2::geom_hline(yintercept = 0, linetype = zero_lty) +
      ggplot2::geom_line(linewidth = lwd_main, color = col_line, na.rm = TRUE) +
      ggplot2::geom_point(size = pt_size, shape = pt_shape, color = col_pts, na.rm = TRUE) +
      ggplot2::geom_text(data = df_na, aes(time, y = 0, label = "NA"),
                         inherit.aes = FALSE, size = 3, vjust = -0.4) +
      ggplot2::labs(x = "Time", y = ttl$ylab,
                    title = ttl$title, subtitle = ttl$subtitle) +
      ggplot2::scale_x_continuous(breaks = brk5) +
      thm +
      ggplot2::theme(plot.title = ggplot2::element_text(color = ttl_col))
    if (add.loess) {
      g <- g + ggplot2::geom_smooth(method = "loess", se = FALSE,
                                    span = span, linewidth = 0.7, color = "black")
    }
    g
  }

  p_resB <- mk_res_panel(
    d_resB, d_naB,
    ttl = list(ylab = "B residuals",
               title = paste0("Bias p-val: ", biasB.p),
               subtitle = NULL),
    ttl_col = colmainB
  )
  p_resF <- mk_res_panel(
    d_resF, d_naF,
    ttl = list(ylab = "F residuals",
               title = paste0("Bias p-val: ", biasF.p),
               subtitle = NULL),
    ttl_col = colmainF
  )

  # ===== ACF panels (osar.acf.plot behavior + dynamic label position) ========
  col_acf_ci <- "blue"  # dashed CI lines
  mk_acf_panel <- function(dacf, pval, ylab_txt, ttl_col) {
    ci <- 1.96 / sqrt(unique(dacf$n))
    # significant lags (exclude lag 0)
    lags_pos <- dacf$lag[dacf$lag > 0]
    acf_pos  <- dacf$acf[dacf$lag > 0]
    sig_idx  <- which(abs(acf_pos) > ci)
    sig_txt  <- if (length(sig_idx) > 0) paste0("lag.signf: ", paste0(lags_pos[sig_idx], collapse = ",")) else ""

    g <- ggplot2::ggplot(dacf, aes(x = lag, y = acf)) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_hline(yintercept = c(-ci, ci), linetype = 2, color = col_acf_ci) +
      ggplot2::geom_segment(aes(xend = lag, y = 0, yend = acf), linewidth = 0.7) +
      ggplot2::labs(x = "Lag", y = ylab_txt,
                    title = paste0("LBox p-val: ", pval)) +
      thm +
      ggplot2::theme(plot.title = ggplot2::element_text(color = ttl_col))

    # dynamic, safe top-right label (never on borders)
    if (nzchar(sig_txt)) {
      xr <- range(dacf$lag, na.rm = TRUE)
      yr <- range(c(dacf$acf, -ci, ci), na.rm = TRUE)
      dx <- xr[2] - xr[1]; dy <- yr[2] - yr[1]
      x_pad <- max(0.03 * dx, 0.15)   # at least 0.15 lag units
      y_pad <- max(0.06 * dy, 0.05)   # at least 0.05 ACF units
      g <- g + ggplot2::annotate("text",
                                 x = xr[2] - x_pad, y = yr[2] - y_pad,
                                 label = sig_txt, hjust = 1, vjust = 1,
                                 color = "red", size = 3)
    }
    g
  }
  p_acfB <- mk_acf_panel(d_acfB, LBoxB.p, "B ACF", ttl_col = colLBoxB)
  p_acfF <- mk_acf_panel(d_acfF, LBoxF.p, "F ACF", ttl_col = colLBoxF)

  # ===== QQ panels (fitted μ+σz line only; title color by p-value) ============
  mk_qq_panel <- function(dqq, pval, ttl_col) {
    ggplot2::ggplot(dqq, aes(theoretical, sample)) +
      ggplot2::geom_line(aes(x = ref_x, y = ref_y), linewidth = 0.7) +
      ggplot2::geom_point(size = 1.8) +
      ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles",
                    title = paste0("Shapiro p-val: ", pval)) +
      thm +
      ggplot2::theme(plot.title = ggplot2::element_text(color = ttl_col))
  }
  p_qqB <- mk_qq_panel(d_qqB, shapiroB.p, ttl_col = colShapB)
  p_qqF <- mk_qq_panel(d_qqF, shapiroF.p, ttl_col = colShapF)

  # ===== Compose grid with vertical row gaps =================================
  # Use spacer rows + relative heights to add vertical space between rows
  row_gap_rel <- 0.14  # adjust this up/down to change the gap between rows

  if (plot.data) {
    grid <- (p_logB + p_logF) /
      patchwork::plot_spacer() /
      (p_resB + p_resF) /
      patchwork::plot_spacer() /
      (p_acfB + p_acfF) /
      patchwork::plot_spacer() /
      (p_qqB  + p_qqF)

    grid <- grid + patchwork::plot_layout(heights = c(1, row_gap_rel, 1, row_gap_rel, 1, row_gap_rel, 1))
  } else {
    grid <- (p_resB + p_resF) /
      patchwork::plot_spacer() /
      (p_acfB + p_acfF) /
      patchwork::plot_spacer() /
      (p_qqB  + p_qqF)

    grid <- grid + patchwork::plot_layout(heights = c(1, row_gap_rel, 1, row_gap_rel, 1))
  }

  grid <- grid +
    patchwork::plot_annotation(caption = stamp) &
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 1, size = 8))

  return(grid)
}
